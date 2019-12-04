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
  if (nrow(A) < 3) {
    tsmessage("Graph too small, using random initialization instead")
    return(rand_init(nrow(A), ndim))
  }

  tsmessage("Initializing from Laplacian Eigenmap")
  # Equivalent to: D <- diag(colSums(A)); M <- solve(D) %*% A
  # This effectively row-normalizes A: colSums is normally faster than rowSums
  # and because A is symmetric, they're equivalent
  M <- A / colSums(A)
  connected <- connected_components(M)
  if (connected$n_components > 1) {
    tsmessage(
      "Found ", connected$n_components, " connected components, ",
      "initializing each component separately"
    )
    fn_name <- as.character(match.call()[[1]])
    return(subgraph_init(fn_name, connected,
      A = A, ndim = ndim,
      verbose = verbose
    ))
  }

  res <- NULL
  k <- ndim + 1
  n <- nrow(M)
  suppressWarnings(
    res <- tryCatch(RSpectra::eigs(M,
      k = k, which = "LM",
      opt = list(tol = 1e-4)
    ),
    error = function(c) {
      NULL
    }
    )
  )

  if (is.null(res) || ncol(res$vectors) < ndim) {
    message(
      "Laplacian Eigenmap failed to converge, ",
      "using random initialization instead"
    )
    return(rand_init(n, ndim))
  }

  vecs <- as.matrix(res$vectors[, 2:(ndim + 1)])
  Re(vecs)
}

# Use a normalized Laplacian.
normalized_laplacian_init <- function(A, ndim = 2, verbose = FALSE) {
  if (nrow(A) < 3) {
    tsmessage("Graph too small, using random initialization instead")
    return(rand_init(nrow(A), ndim))
  }
  tsmessage("Initializing from normalized Laplacian")
  connected <- connected_components(A)
  if (connected$n_components > 1) {
    tsmessage(
      "Found ", connected$n_components, " connected components, ",
      "initializing each component separately"
    )
    fn_name <- as.character(match.call()[[1]])
    return(subgraph_init(fn_name, connected,
      A = A, ndim = ndim,
      verbose = verbose
    ))
  }

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
  opt <- list(tol = 1e-4)
  suppressWarnings(
    res <- tryCatch(RSpectra::eigs_sym(L, k = k, which = "SM", opt = opt),
      error = function(c) {
        NULL
      }
    )
  )
  if (is.null(res) || ncol(res$vectors) < ndim) {
    suppressWarnings(
      res <- tryCatch(RSpectra::eigs_sym(L,
        k = k, which = "LM", sigma = 0,
        opt = opt
      ),
      error = function(c) {
        NULL
      }
      )
    )
    if (is.null(res) || ncol(res$vectors) < ndim) {
      message(
        "Spectral initialization failed to converge, ",
        "using random initialization instead"
      )
      return(rand_init(n, ndim))
    }
  }
  vec_indices <- rev(order(res$values, decreasing = TRUE)[1:ndim])
  as.matrix(Re(res$vectors[, vec_indices]))
}

# Use irlba's partial_eigen instead of RSpectra
irlba_normalized_laplacian_init <- function(A, ndim = 2, verbose = FALSE) {
  if (nrow(A) < 3) {
    tsmessage("Graph too small, using random initialization instead")
    return(rand_init(nrow(A), ndim))
  }
  tsmessage("Initializing from normalized Laplacian (using irlba)")

  n <- nrow(A)
  Dsq <- sqrt(Matrix::colSums(A))
  L <- -Matrix::t(A / Dsq) / Dsq
  Matrix::diag(L) <- 1 + Matrix::diag(L)

  k <- ndim + 1

  suppressWarnings(
    res <- tryCatch(res <- irlba::partial_eigen(L,
      n = k, symmetric = TRUE,
      smallest = TRUE, tol = 1e-3,
      maxit = 1000, verbose = TRUE
    ),
    error = function(c) {
      NULL
    }
    )
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
  if (nrow(A) < 3) {
    tsmessage("Graph too small, using random initialization instead")
    return(rand_init(nrow(A), ndim))
  }
  tsmessage("Initializing from normalized Laplacian + noise")
  connected <- connected_components(A)
  if (connected$n_components > 1) {
    tsmessage(
      "Found ", connected$n_components, " connected components, ",
      "initializing each component separately"
    )
    fn_name <- as.character(match.call()[[1]])
    return(subgraph_init(fn_name, connected,
      A = A, ndim = ndim,
      verbose = verbose
    ))
  }
  coords <- normalized_laplacian_init(A, ndim, verbose = FALSE)
  expansion <- 10.0 / max(abs(coords))
  (coords * expansion) + matrix(stats::rnorm(n = prod(dim(coords)), sd = 0.0001),
    ncol = ndim
  )
}

irlba_spectral_init <- function(A, ndim = 2, verbose = FALSE) {
  if (nrow(A) < 3) {
    tsmessage("Graph too small, using random initialization instead")
    return(rand_init(nrow(A), ndim))
  }
  tsmessage("Initializing from normalized Laplacian (using irlba) + noise")

  coords <- irlba_normalized_laplacian_init(A, ndim, verbose = FALSE)
  expansion <- 10.0 / max(coords)
  (coords * expansion) + matrix(stats::rnorm(n = prod(dim(coords)), sd = 0.001),
    ncol = ndim
  )
}

# Recursively calls the spectral initialization function named fn_name
# for each subgraph specified by connected
subgraph_init <- function(fn_name, connected, A, ndim = 2, verbose = FALSE) {
  init <- NULL
  for (i in 1:connected$n_components) {
    subg_idx <- connected$labels == i - 1
    subg <- A[subg_idx, subg_idx]
    tsmessage("Initializing subcomponent of size ", nrow(subg))
    init_conn <- do.call(fn_name, list(
      A = subg, ndim = ndim,
      verbose = verbose
    ))
    if (is.null(init)) {
      init <- init_conn
    }
    else {
      init <- rbind(init, init_conn)
    }
  }
  init
}

# Return the number of connected components in a graph (respresented as a
# sparse matrix).
connected_components <- function(X) {
  Xt <- Matrix::t(X)
  connected_components_undirected(nrow(X), Xt@i, Xt@p, X@i, X@p)
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

# Rescale embedding so that the standard deviation is the specified value.
# Default gives initialization like t-SNE, but not random. Large initial
# distances lead to small gradients, and hence small updates, so should be
# avoided
shrink_coords <- function(X, sdev = 1e-4) {
  scale(X, scale = apply(X, 2, stats::sd) / sdev)
}

# PCA
pca_init <- function(X, ndim = 2, center = TRUE, verbose = FALSE) {
  tsmessage("Initializing from PCA")
  pca_scores(X, ncol = ndim, center = center, verbose = verbose)
}


# Calculates a matrix containing the first ncol columns of the PCA scores.
# Returns the score matrix unless ret_extra is TRUE, in which case a list
# is returned also containing the eigenvalues
pca_scores <- function(X, ncol = min(dim(X)), center = TRUE, ret_extra = FALSE,
                       verbose = FALSE) {
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
    return(scores)
  }

  # irlba warns about using too large a percentage of total singular value
  # so don't use if dataset is small compared to ncol
  if (ncol < 0.5 * min(dim(X))) {
    return(irlba_scores(X,
      ncol = ncol, center = center, ret_extra = ret_extra,
      verbose = verbose
    ))
  }

  svd_scores(X = X, ncol = ncol, center = center, ret_extra = ret_extra, verbose = verbose)
}

# Get scores by SVD
svd_scores <- function(X, ncol = min(dim(X)), center = TRUE, ret_extra = FALSE,
                       verbose = FALSE) {
  # need extra data if we want to re-apply PCA to new points in umap_transform
  rotation <- NULL
  xcenter <- NULL

  X <- scale(X, center = center, scale = FALSE)
  # do SVD on X directly rather than forming covariance matrix
  s <- svd(X, nu = ncol, nv = ifelse(ret_extra, ncol, 0))
  D <- diag(c(s$d[1:ncol]), ncol, ncol)
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
  if (ret_extra) {
    rotation <- s$v
    xcenter <- attr(X, "scaled:center")
  }

  if (ret_extra) {
    list(
      scores = scores,
      lambda = lambda[1:ncol],
      rotation = rotation,
      center = xcenter
    )
  }
  else {
    scores
  }
}

# Get PCA scores via irlba
irlba_scores <- function(X, ncol, center = TRUE, ret_extra = FALSE, verbose = FALSE) {
  res <- irlba::prcomp_irlba(X,
    n = ncol, retx = TRUE, center = center,
    scale = FALSE
  )
  if (verbose) {
    varex <- sum(res$sdev[1:ncol]^2) / res$totalvar
    tsmessage(
      "PCA: ", ncol, " components explained ", formatC(varex * 100),
      "% variance"
    )
  }
  if (ret_extra) {
    list(scores = res$x, rotation = res$rotation, center = res$center)
  }
  else {
    res$x
  }
}

init_is_spectral <- function(init) {
  res <- pmatch(tolower(init), c(
    "normlaplacian", "spectral", "laplacian",
    "inormlaplacian", "ispectral"
  ))
  length(res) > 0 && !is.na(res)
}

rand_nbr_graph <- function(n_vertices, n_nbrs, val) {
  nn_to_sparse(rand_nbr_idx(n_vertices, n_nbrs),
    val = val,
    max_nbr_id = n_vertices
  )
}

rand_nbr_idx <- function(n_vertices, n_nbrs) {
  idx <- matrix(nrow = n_vertices, ncol = n_nbrs)
  nv1 <- n_vertices - 1
  for (i in 1:n_vertices) {
    ids <- sample.int(nv1, n_nbrs)
    id_sel <- ids >= 1
    ids[id_sel] <- ids[id_sel] + 1
    idx[i, ] <- ids
  }
  idx
}

# V: the current affinity graph
# n_pos: number of neighbors to retain per item
# n_neg: number of "negative" (i.e. non-)neighbors per item
# pos_affinity: value for the positive affinity (associated with nbrs)
# neg_affinity: value for the negative affinity (associated with neg nbrs)
approx_affinity_graph <- function(V, n_neg,
                                  pos_affinity = 1, neg_affinity = 0.1,
                                  verbose = FALSE) {
  pos_V <- V
  pos_V@x <- rep(pos_affinity, length(pos_V@x))
  pos_V <- 0.5 * (pos_V + Matrix::t(pos_V))

  neg_V <- rand_nbr_graph(nrow(pos_V), n_nbrs = n_neg, val = neg_affinity)
  neg_V <- 0.5 * (neg_V + Matrix::t(neg_V))

  # the cleanup below will ensure that where the same value got a pos and neg
  # affinity it will end up positive
  graph <- pos_V + neg_V

  # clamp small values to neg_affinity
  graph@x[graph@x < pos_affinity] <- neg_affinity
  # and large values to pos_affinity
  graph@x <- pmin(graph@x, pos_affinity)

  Matrix::drop0(graph)
}


# Initialize using a spectral decomposition of an "approximate global" graph
# Uses the same graph as standard UMAP, but with each entry set to 1. A measure
# of global structure is added by randomly setting some of the remaining zero
# to a smaller value (0.1 in this case).
# This routine is inspired by some ideas in
# 2-D Embedding of Large and High-dimensional Data with Minimal Memory and Computational Time Requirements
# Witold Dzwinel, Rafal Wcislo, Stan Matwin
# https://arxiv.org/abs/1902.01108
#
# Randomized Near Neighbor Graphs, Giant Components, and Applications in Data Science
# George C. Linderman, Gal Mishne, Yuval Kluger, Stefan Steinerberger
# https://arxiv.org/abs/1711.04712
agspectral_init <- function(V, n_neg_nbrs, pos_affinity = 1, neg_affinity = 0.1,
                            ndim = 2, verbose = FALSE) {
  graph <- approx_affinity_graph(V, n_neg_nbrs,
    pos_affinity = pos_affinity,
    neg_affinity = neg_affinity,
    verbose = verbose
  )
  spectral_init(graph, ndim = ndim, verbose = verbose)
}
