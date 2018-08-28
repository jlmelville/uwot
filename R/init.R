# Laplacian Eigenmap (Belkin & Niyogi, 2002)
# Original formulation solves the generalized eigenvalue problem of the
# unnormalized graph Laplacian: L v = lambda D v, where L = D - A
# and uses the bottom eigenvectors v that result
# (ignoring the constant eigenvector associated with the smallest eigenvalue).
#
# This is equivalent to using the top eigenvectors from the usual
# eigendecomposition of a row-normalized Laplacian P = D^-1 A: P v = lambda' v
# so we don't need to depend on an external package for generalized eigenvalues.
# Note that while the eigenvectors are the same, the eigenvalues are
# different: lambda' = 1 - lambda, but we don't use them with Laplacian
# Eigenmaps anyway.
#
# As we only need to calculate the top ndim + 1 eigenvectors (i.e. normally 3)
# it's incredibly wasteful to calculate all of them.
# A must be symmetric and positive semi definite, but not necessarily
# normalized in any specific way.
#' @import Matrix
laplacian_eigenmap <- function(A, ndim = 2, verbose = FALSE) {
  tsmessage("Initializing from Laplacian Eigenmap")
  # Equivalent to: D <- diag(colSums(A)); M <- solve(D) %*% A
  # This effectively row-normalizes A: colSums is normally faster than rowSums
  # and because A is symmetric, they're equivalent
  M <- A / colSums(A)

  eig_res <- tryCatch(RSpectra::eigs(M, k = ndim + 1),
    error = function(c) {
      NULL
    }
  )

  if (is.null(eig_res) || ncol(eig_res$vectors) < ndim + 1) {
    message(
      "Laplacian Eigenmap failed to converge, ",
      "using random initialization instead"
    )
    return(rand_init(nrow(A), ndim))
  }
  vecs <- as.matrix(eig_res$vectors[, 2:(ndim + 1)])
  Re(vecs)
}

# Use a normalized Laplacian.
normalized_laplacian_init <- function(A, ndim = 2, verbose = FALSE) {
  tsmessage("Initializing from normalized Laplacian")

  n <- nrow(A)
  # Normalized Laplacian: clear and close to UMAP code, but very slow in R
  # I <- diag(1, nrow = n, ncol = n)
  # D <- diag(1 / sqrt(colSums(A)))
  # L <- I - D %*% A %*% D

  # A lot faster (order of magnitude when n = 1000)
  Dsq <- sqrt(Matrix::colSums(A))
  L <- -Matrix::t(A / Dsq) / Dsq
  Matrix::diag(L) <- 1 + Matrix::diag(L)

  k <- ndim + 1
  ncv <- max(2 * k + 1, floor(sqrt(n)))
  opt <- list(
    ncv = ncv,
    maxitr = 5 * n,
    tol = 1e-4
  )
  res <- tryCatch(RSpectra::eigs_sym(L, k = k, which = "SM", opt = opt),
    error = function(c) {
      NULL
    }
  )
  if (is.null(res) || ncol(res$vectors) < ndim) {
    message(
      "Spectral initialization failed to converge, ",
      "using random initialization instead"
    )
    return(rand_init(n, ndim))
  }
  vec_indices <- rev(order(res$values, decreasing = TRUE)[1:ndim])
  as.matrix(Re(res$vectors[, vec_indices]))
}

# Default UMAP initialization
# spectral decomposition of the normalized Laplacian + some noise
spectral_init <- function(A, ndim = 2, verbose = FALSE) {
  tsmessage("Initializing from normalized Laplacian + noise")
  coords <- normalized_laplacian_init(A, ndim, verbose = FALSE)
  expansion <- 10.0 / max(coords)
  (coords * expansion) + matrix(stats::rnorm(n = prod(dim(coords)), sd = 0.001),
    ncol = ndim
  )
}

# UMAP random initialization: uniform between +10 and -10 along each axis
rand_init <- function(n, ndim, verbose = FALSE) {
  tsmessage("Initializing from uniform random")
  matrix(stats::runif(n = n * ndim, min = -10, max = 10), ncol = ndim)
}

# LargeVis random initialization: Gaussian with sd 1e-4 (like t-SNE)
rand_init_lv <- function(n, ndim, verbose = FALSE) {
  tsmessage("Initializing from random Gaussian with sd = 1e-4")
  matrix(stats::rnorm(ndim * n, sd = 1e-4), n)
}

# PCA but then scale the vectors to a t-SNE-like stdev of 1e-4
scaled_pca <- function(X, ndim = 2, verbose = FALSE) {
  tsmessage("Initializing from scaled PCA")
  scores <- pca_scores(X, ncol = ndim, verbose = verbose)
  scale(scores, scale = apply(scores, 2, stats::sd) / 1e-4)
}

# PCA
pca_init <- function(X, ndim = 2, verbose = FALSE) {
  tsmessage("Initializing from PCA")
  pca_scores(X, ncol = ndim, verbose = verbose)
}


# Calculates a matrix containing the first ncol columns of the PCA scores.
# Returns the score matrix unless ret_extra is TRUE, in which case a list
# is returned also containing the eigenvalues
pca_scores <- function(X, ncol = min(dim(X)), ret_extra = FALSE,
                       verbose = FALSE) {
  # irlba warns about using too large a percentage of total singular value
  # so don't use if dataset is small compared to ncol
  if (ncol < 0.5 * min(dim(X))) {
    return(irlba_scores(X, ncol = ncol))
  }

  if (methods::is(X, "dist")) {
    res_mds <- stats::cmdscale(X, x.ret = TRUE, eig = TRUE, k = ncol)

    if (ret_extra || verbose) {
      lambda <- res_mds$eig
      varex <- sum(lambda[1:ncol]) / sum(lambda)
      tsmessage(
        "PCA (using classical MDS): ", ncol, " components explained ",
        formatC(varex * 100), "% variance"
      )
    }
    scores <- res_mds$points
  }
  else {
    X <- scale(X, center = TRUE, scale = FALSE)
    # do SVD on X directly rather than forming covariance matrix
    s <- svd(X, nu = ncol, nv = 0)
    D <- diag(c(s$d[1:ncol]))
    if (verbose || ret_extra) {
      # calculate eigenvalues of covariance matrix from singular values
      lambda <- (s$d^2) / (nrow(X) - 1)
      varex <- sum(lambda[1:ncol]) / sum(lambda)
      tsmessage(
        "PCA: ", ncol, " components explained ", formatC(varex * 100),
        "% variance"
      )
    }
    scores <- s$u %*% D
  }

  if (ret_extra) {
    list(
      scores = scores,
      lambda = lambda[1:ncol]
    )
  }
  else {
    scores
  }
}

# Get PCA scores via irlba
irlba_scores <- function(X, ncol) {
  irlba::prcomp_irlba(X, n = ncol, retx = TRUE, center = TRUE, scale = FALSE)$x
}
