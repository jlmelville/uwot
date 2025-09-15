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
laplacian_eigenmap <- function(A, ndim = 2, verbose = FALSE, force_irlba = FALSE) {
  if (rspectra_is_installed() && !force_irlba) {
    coords <- rspectra_laplacian_eigenmap(A, ndim, verbose = verbose)
  } else {
    coords <- irlba_laplacian_eigenmap(A, ndim, verbose = verbose)
  }
  coords
}

rspectra_laplacian_eigenmap <- function(A, ndim = 2, verbose = FALSE) {
  if (nrow(A) < 3) {
    tsmessage("Graph too small, using random initialization instead")
    return(rand_init(nrow(A), ndim))
  }

  tsmessage("Initializing from Laplacian Eigenmap (via RSpectra)")
  # Equivalent to: D <- diag(colSums(A)); M <- solve(D) %*% A
  # This effectively row-normalizes A: colSums is normally faster than rowSums
  # and because A is symmetric, they're equivalent
  M <- A / colSums(A)
  res <- rspectra_eigs_asym(M, ndim)
  if (is.null(res) ||
    !is.list(res) ||
    !"vectors" %in% names(res) ||
    is.null(res$vectors) ||
    tryCatch(
      is.na(ncol(res$vectors)),
      error = function(e) {
        TRUE
      }
    ) ||
    ncol(res$vectors) < ndim) {
    message(
      "Laplacian Eigenmap failed to converge, ",
      "using random initialization instead"
    )
    n <- nrow(M)
    return(rand_init(n, ndim))
  }

  # return the smallest eigenvalues
  as.matrix(Re(res$vectors[, 2:(ndim + 1)]))
}

irlba_laplacian_eigenmap <- function(A, ndim = 2, verbose = FALSE) {
  if (nrow(A) < 3) {
    tsmessage("Graph too small, using random initialization instead")
    return(rand_init(nrow(A), ndim))
  }
  tsmessage("Initializing from Laplacian Eigenmap (via irlba)")

  lapA <- form_modified_laplacian(A, ret_d = TRUE)
  res <- irlba_spectral_tsvd(lapA$L, ndim + 1)
  if (is.null(res) || ncol(res$vectors) < ndim || !res$converged) {
    message(
      "Laplacian Eigenmap failed to converge, ",
      "using random initialization instead"
    )
    return(rand_init(nrow(A), ndim))
  }
  res <- lapA$Disqrt * res$vectors[, 2:(ndim + 1), drop = FALSE]
  # re-scale the vectors to length 1
  sweep(res, 2, sqrt(colSums(res * res)), `/`)
}

form_normalized_laplacian <- function(A) {
  # Normalized Laplacian: clear and close to UMAP code, but very slow in R
  # I <- diag(1, nrow = n, ncol = n)
  # D <- diag(1 / sqrt(colSums(A)))
  # L <- I - D %*% A %*% D

  # A lot faster (order of magnitude when n = 1000)
  Dsq <- sqrt(Matrix::colSums(A))
  L <- -Matrix::t(A / Dsq) / Dsq
  Matrix::diag(L) <- 1 + Matrix::diag(L)
  L
}

# The symmetrized graph Laplacian (Lsym) but shifted so that:
# the bottom eigenvectors of Lsym correspond to the top singular vectors of
# this matrix (hence can be used with truncated SVD), and the eigenvalues
# are all positive, so we don't lose sign and hence correct eigenvector ordering
# when using the singular values (lambda = 2 - d)
# effectively we form 2I - Lsym = D^-1/2 W D^-1/2 + I
form_modified_laplacian <- function(A, ret_d = FALSE) {
  Dsq <- sqrt(Matrix::colSums(A))
  L <- Matrix::t(A / Dsq) / Dsq
  Matrix::diag(L) <- 1 + Matrix::diag(L)
  if (ret_d) {
    list(L = L, Disqrt = 1 / Dsq)
  } else {
    L
  }
}

# Return the ndim eigenvectors associated with the ndim largest eigenvalues
sort_eigenvectors <- function(eig_res, ndim, decreasing = TRUE) {
  vec_indices <- rev(order(eig_res$values, decreasing = decreasing)[1:ndim])
  as.matrix(Re(eig_res$vectors[, vec_indices]))
}

normalized_laplacian_init <- function(A, ndim = 2, verbose = FALSE, force_irlba = FALSE) {
  if (rspectra_is_installed() && !force_irlba) {
    coords <- rspectra_normalized_laplacian_init(A, ndim, verbose = verbose)
  } else {
    coords <- irlba_normalized_laplacian_init(A, ndim, verbose = verbose)
  }
  coords
}

rspectra_normalized_laplacian_init <- function(A, ndim = 2, verbose = FALSE) {
  if (nrow(A) < 3) {
    tsmessage("Graph too small, using random initialization instead")
    return(rand_init(nrow(A), ndim))
  }
  tsmessage("Initializing from normalized Laplacian")

  L <- form_normalized_laplacian(A)
  res <- rspectra_eigs_sym(L, ndim, verbose = verbose)

  if (is.null(res) ||
      !is.list(res) ||
      !"vectors" %in% names(res) ||
      is.null(res$vectors) ||
      tryCatch(
        is.na(ncol(res$vectors)),
        error = function(e) {
          TRUE
        }
      ) ||
      ncol(res$vectors) < ndim) {
    message(
      "Spectral initialization failed to converge, ",
      "using random initialization instead"
    )
    n <- nrow(A)
    return(rand_init(n, ndim))
  }
  sort_eigenvectors(res, ndim)
}

# maybe this should become an option one day
# reminder to me for the next time I experiment: Shift-invert just doesn't work
# in all cases: causes MNIST (k = 15) to hang forever even with a guess for
# initvec of D^1/2, no matter which shift value (or maxitr or tol value) is
# used. But we can form the shifted Laplacian like we do for tsvd approaches
# and look for the LM eigenvalues (plus initialize with D^1/2 as a guess)
# However the tolerance needs to be lower for similar quality of output
# and most initializations are < 1 second anyway, so we don't gain much from
# speeding them up. For slower initializations (e.g. `tomoradar`, `mammoth`)
# there is a small speed up, but not down to a few seconds
rspectra_normalized_laplacian_init_shift_inv <- function(A, ndim = 2, verbose = FALSE) {
  if (nrow(A) < 3) {
    tsmessage("Graph too small, using random initialization instead")
    return(rand_init(nrow(A), ndim))
  }
  tsmessage("Initializing from normalized Laplacian")

  L <- form_modified_laplacian(A)
  initvec <- sqrt(Matrix::colSums(A))
  res <- rspectra_eigs_shift_sym(L,
    ndim,
    verbose = verbose,
    initvec = initvec,
    tol = 1e-6
  )
  norm_initvec <- initvec / sqrt(sum(initvec * initvec))
  if (is.null(res) ||
    !is.list(res) ||
    !"vectors" %in% names(res) ||
    is.null(res$vectors) ||
    tryCatch(
      is.na(ncol(res$vectors)),
      error = function(e) {
        TRUE
      }
    ) ||
    ncol(res$vectors) < ndim) {
    message(
      "Spectral initialization failed to converge, ",
      "using random initialization instead"
    )
    n <- nrow(A)
    return(rand_init(n, ndim))
  }
  sort_eigenvectors(res, ndim, decreasing = FALSE)
}

# Use a normalized Laplacian and use truncated SVD
irlba_tsvd_normalized_laplacian_init <- function(A, ndim = 2, verbose = FALSE) {
  if (nrow(A) < 3) {
    tsmessage("Graph too small, using random initialization instead")
    return(rand_init(nrow(A), ndim))
  }
  tsmessage("Initializing from normalized Laplacian")

  L <- form_modified_laplacian(A)
  res <- irlba_spectral_tsvd(L, ndim + 1)
  if (is.null(res) || ncol(res$vectors) < ndim || !res$converged) {
    message(
      "Spectral initialization failed to converge, ",
      "using random initialization instead"
    )
    n <- nrow(A)
    return(rand_init(n, ndim))
  }
  res$vectors[, 2:(ndim + 1), drop = FALSE]
}

irlba_spectral_tsvd <- function(L, n, iters = 1000) {
  irlba_args <- list(
    A = L,
    nv = n,
    nu = 0,
    maxit = iters
  )
  suppressWarnings(res <- tryCatch(
    do.call(irlba::irlba, irlba_args),
    error = function(c) {
      irlba_args$fastpath <- FALSE
      do.call(irlba::irlba, irlba_args)
    }
  ))
  list(
    vectors = res$v,
    values = 2.0 - res$d,
    converged = res$iter != iters
  )
}

irlba_eigs_asym <- function(L, ndim) {
  irlba_args <- list(
    x = L,
    n = ndim + 1,
    symmetric = FALSE,
    smallest = TRUE,
    tol = 1e-3,
    maxit = 1000
  )
  suppressWarnings(res <- tryCatch(
    do.call(irlba::partial_eigen, irlba_args),
    error = function(e) {
      irlba_args$fastpath <- FALSE
      tryCatch(
        do.call(irlba::partial_eigen, irlba_args),
        error = function(e) {
          NULL
        }
      )
    }
  ))
  if (!is.null(res)) {
    res$values <- sqrt(res$values)
  }
  res
}

irlba_eigs_sym <- function(L, ndim, smallest = TRUE) {
  irlba_args <- list(
    x = L,
    n = ndim + 1,
    symmetric = TRUE,
    smallest = smallest,
    tol = 1e-3,
    maxit = 1000
  )
  suppressWarnings(res <- tryCatch(
    do.call(irlba::partial_eigen, irlba_args),
    error = function(e) {
      irlba_args$fastpath <- FALSE
      tryCatch(
        do.call(irlba::partial_eigen, irlba_args),
        error = function(e) {
          NULL
        }
      )
    }
  ))
  res
}

# Use irlba's partial_eigen instead of RSpectra
irlba_normalized_laplacian_init <- function(A, ndim = 2, verbose = FALSE) {
  if (nrow(A) < 3) {
    tsmessage("Graph too small, using random initialization instead")
    return(rand_init(nrow(A), ndim))
  }
  tsmessage("Initializing from normalized Laplacian (using irlba)")

  # Using the normalized Laplacian and looking for smallest eigenvalues does
  # not work well with irlba's partial_eigen routine, so form the shifted
  # Laplacian and look for largest eigenvalues
  L <- form_modified_laplacian(A)
  res <- irlba_eigs_sym(L, ndim, smallest = FALSE)
  if (is.null(res) ||
    !is.list(res) ||
    !"vectors" %in% names(res) ||
    is.null(res$vectors) ||
    tryCatch(
      is.na(ncol(res$vectors)),
      error = function(e) {
        TRUE
      }
    ) ||
    ncol(res$vectors) < ndim) {
    message(
      "Spectral initialization failed to converge, ",
      "using random initialization instead"
    )
    n <- nrow(A)
    return(rand_init(n, ndim))
  }
  # shift back the eigenvalues
  res$values <- 2.0 - res$values
  sort_eigenvectors(res, ndim)
}


# Default UMAP initialization
# spectral decomposition of the normalized Laplacian + some noise
spectral_init <- function(A, ndim = 2, verbose = FALSE, force_irlba = FALSE) {
  if (nrow(A) < 3) {
    tsmessage("Graph too small, using random initialization instead")
    return(rand_init(nrow(A), ndim))
  }
  if (rspectra_is_installed() && !force_irlba) {
    tsmessage("Initializing from normalized Laplacian + noise (using RSpectra)")
    coords <- rspectra_normalized_laplacian_init(A, ndim, verbose = FALSE)
  } else {
    tsmessage("Initializing from normalized Laplacian + noise (using irlba)")
    coords <- irlba_tsvd_normalized_laplacian_init(A, ndim, verbose = FALSE)
  }
  scale_and_jitter(coords, max_coord = 10.0, sd = 0.0001)
}

irlba_spectral_init <- function(A, ndim = 2, verbose = FALSE) {
  if (nrow(A) < 3) {
    tsmessage("Graph too small, using random initialization instead")
    return(rand_init(nrow(A), ndim))
  }
  tsmessage("Initializing from normalized Laplacian (using irlba) + noise")

  coords <- irlba_normalized_laplacian_init(A, ndim, verbose = FALSE)
  scale_and_jitter(coords, max_coord = 10.0, sd = 0.0001)
}

# Scales coords so that the largest absolute coordinate is 10.0 then jitters by
# adding gaussian noise with mean 0 and standard deviation sd
scale_and_jitter <- function(coords, max_coord = 10.0, sd = 0.0001) {
  expansion <- 10.0 / max(abs(coords))
  (coords * expansion) + matrix(stats::rnorm(n = prod(dim(coords)), sd = sd),
    ncol = ncol(coords)
  )
}

# Return the number of connected components in a graph (represented as a
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
# avoided.
scale_coords <- function(X, sdev = 1e-4, verbose = FALSE) {
  if (is.character(sdev) && sdev == "range") {
    # #99: range scale coordinates like python UMAP does
    tsmessage("Range-scaling initial input columns to 0-10")
    return(apply(X, 2, range_scale, max = 10.0))
  }

  if (is.null(sdev)) {
    return(X)
  }
  tsmessage("Scaling init to sdev = ", sdev)
  scale_factor <- apply(X, 2, stats::sd)
  scale(X, scale = scale_factor / sdev)
}

# PCA
# Calculates a matrix containing the first ndim columns of the PCA scores.
# Returns the score matrix unless ret_extra is TRUE, in which case a list
# is returned also containing the eigenvalues
pca_init <- function(X, ndim = min(dim(X)), center = TRUE, ret_extra = FALSE,
                     pca_method = "auto", verbose = FALSE) {
  if (inherits(X, "dist")) {
    res_mds <- stats::cmdscale(X, x.ret = TRUE, eig = TRUE, k = ndim)
    if (ret_extra || verbose) {
      lambda <- res_mds$eig
      lambda[lambda < 0] <- 0
      varex <- sum(lambda[1:ndim]) / sum(lambda)
      tsmessage(
        "PCA (using classical MDS): ", ndim, " components explained ",
        formatC(varex * 100), "% variance"
      )
    }
    scores <- res_mds$points
    return(scores)
  }

  # irlba warns about using too large a percentage of total singular value
  # so don't use if dataset is small compared to ndim
  if (pca_method == "auto") {
    if (ndim < 0.5 * min(dim(X))) {
      pca_method <- "irlba"
    } else {
      pca_method <- "svd"
    }
  }

  if (pca_method == "bigstatsr") {
    if (!bigstatsr_is_installed()) {
      warning(
        "PCA via bigstatsr requires the 'bigstatsr' package. ",
        "Please install it. Falling back to 'irlba'"
      )
      pca_method <- "irlba"
    }
  }

  tsmessage("Using '", pca_method, "' for PCA")
  pca_fun <- switch(pca_method,
    irlba = irlba_scores,
    svdr = irlba_svdr_scores,
    svd = svd_scores,
    bigstatsr = bigstatsr_scores,
    stop("BUG: unknown svd method '", pca_method, "'")
  )

  do.call(pca_fun, list(
    X = X,
    ncol = ndim,
    center = center,
    ret_extra = ret_extra,
    verbose = verbose
  ))
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
  } else {
    scores
  }
}

# Get PCA scores via irlba
irlba_scores <- function(X, ncol, center = TRUE, ret_extra = FALSE, verbose = FALSE) {
  if (is.logical(X)) {
    tsmessage("Converting logical input to numeric for PCA initialization")
    # convert logical matrix to numeric
    X <- X * 1
  }
  res <- irlba::prcomp_irlba(X,
    n = ncol, retx = TRUE, center = center, scale. = FALSE
  )
  report_varex(res, verbose)
  if (ret_extra) {
    list(scores = res$x, rotation = res$rotation, center = res$center)
  } else {
    res$x
  }
}

report_varex <- function(res, verbose = FALSE) {
  if (verbose) {
    ncol <- ncol(res$rotation)
    varex <- sum(res$sdev[1:ncol]^2) / res$totalvar
    tsmessage(
      "PCA: ",
      ncol,
      " components explained ",
      formatC(varex * 100),
      "% variance"
    )
  }
}

# This function taken from irlba and modified to use irlba::svdr rather
# than irlba::irlba
prcomp_rsvd <- function(x, n = 3, retx = TRUE, center = TRUE, scale. = FALSE,
                        ...) {
  a <- names(as.list(match.call()))
  ans <- list(scale = scale.)
  if ("tol" %in% a) {
    warning("The `tol` truncation argument from `prcomp` is not supported by\n`prcomp_rsvd`. If specified, `tol` is passed to the `irlba` function to\ncontrol that algorithm's convergence tolerance. See `?prcomp_irlba` for help.")
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  args <- list(x = x, k = n)
  if (is.logical(center)) {
    if (center) {
      args$center <- colMeans(x)
    }
  } else {
    args$center <- center
  }
  if (is.logical(scale.)) {
    if (is.numeric(args$center)) {
      f <- function(i) {
        sqrt(sum((x[, i] - args$center[i])^2) / (nrow(x) -
          1L))
      }
      scale. <- vapply(seq(ncol(x)), f, pi, USE.NAMES = FALSE)
      if (ans$scale) {
        ans$totalvar <- ncol(x)
      } else {
        ans$totalvar <- sum(scale.^2)
      }
    } else {
      if (ans$scale) {
        scale. <- apply(x, 2L, function(v) {
          sqrt(sum(v^2) / max(
            1,
            length(v) - 1L
          ))
        })
        f <- function(i) {
          sqrt(sum((x[, i] / scale.[i])^2) / (nrow(x) -
            1L))
        }
        ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi,
          USE.NAMES = FALSE
        )^2)
      } else {
        f <- function(i) sum(x[, i]^2) / (nrow(x) - 1L)
        ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi,
          USE.NAMES = FALSE
        ))
      }
    }
    if (ans$scale) {
      args$scale <- scale.
    }
  } else {
    args$scale <- scale.
    f <- function(i) {
      sqrt(sum((x[, i] / scale.[i])^2) / (nrow(x) -
        1L))
    }
    ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES = FALSE))
  }
  if (!missing(...)) {
    args <- c(args, list(...))
  }

  s <- do.call(irlba::svdr, args = args)
  ans$sdev <- s$d / sqrt(max(1, nrow(x) - 1))
  ans$rotation <- s$v
  colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)),
    sep = ""
  )
  ans$center <- args$center
  if (retx) {
    ans <- c(ans, list(x = sweep(s$u, 2, s$d, FUN = `*`)))
    colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)),
      sep = ""
    )
  }
  class(ans) <- c("irlba_prcomp", "prcomp")
  ans
}

irlba_svdr_scores <-
  function(X,
           ncol,
           center = TRUE,
           ret_extra = FALSE,
           verbose = FALSE) {
    # 5 iterations is the default for scikit-learn TruncatedSVD
    res <- prcomp_rsvd(
      X,
      n = ncol,
      retx = TRUE,
      center = center,
      scale. = FALSE,
      it = 5
    )
    report_varex(res, verbose)
    if (ret_extra) {
      list(
        scores = res$x,
        rotation = res$rotation,
        center = res$center
      )
    } else {
      res$x
    }
  }

init_is_spectral <- function(init) {
  res <- pmatch(tolower(init), c(
    "normlaplacian", "spectral", "laplacian",
    "inormlaplacian", "ispectral", "agspectral",
    "irlba_spectral", "irlba_laplacian"
  ))
  length(res) > 0 && !is.na(res)
}

rand_nbr_graph <- function(n_vertices, n_nbrs, val) {
  nng_to_sparse(rand_nbr_idx(n_vertices, n_nbrs),
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
