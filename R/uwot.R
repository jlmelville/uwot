#' @useDynLib uwot
#' @importFrom Rcpp sourceCpp
NULL

umap <- function(X, n_neighbors = 15, n_components = 2, n_epochs = NULL,
                 alpha = 1, init = "spectral", spread = 1, min_dist = 0.01,
                 set_op_mix_ratio = 1.0, local_connectivity = 1.0,
                 bandwidth = 1.0, gamma = 1.0,
                 negative_sample_rate = 5.0, a = NULL, b = NULL,
                 nn_method = NULL,
                 verbose = getOption("verbose", TRUE)) {

  if (is.null(a) || is.null(b)) {
    ab_res <- find_ab_params(spread = spread, min_dist = min_dist)
    a <- ab_res[1]
    b <- ab_res[2]
  }

  if (methods::is(X, "dist")) {
    n_vertices <- attr(X, "Size")
  }
  else {
    if (methods::is(X, "data.frame")) {
      indexes <- which(vapply(X, is.numeric, logical(1)))
      if (verbose) {
        message("Found ", length(indexes), " numeric columns")
      }
      if (length(indexes) == 0) {
        stop("No numeric columns found")
      }
      X <- as.matrix(X[, indexes])
    }
    n_vertices <- nrow(X)
  }

  if (is.null(nn_method)) {
    if (n_vertices < 4096) {
      if (verbose) {
        tsmessage("Using FNN for neighbor search")
      }
      nn_method = "fnn"
    }
    else {
      if (verbose) {
        tsmessage("Using ANNOY for neighbor search")
      }
      nn_method = "annoy"
    }
  }

  V <- fuzzy_set_union(smooth_knn_distances(nn = find_nn(X, n_neighbors,
                                                         method = nn_method),
                                            local_connectivity = local_connectivity,
                                            bandwidth = bandwidth,
                                            verbose = verbose
                                            )$P)

  if (methods::is(init, "matrix")) {
    if (nrow(init) != n_vertices || ncol(init) != n_components) {
      stop("init matrix does not match necessary configuration for X")
    }
  }
  else {
    if (init == "spectral") {
      embedding <- spectral_init(V, ndim = n_components, use_RSpectra = TRUE,
                                 verbose = verbose)
    }
    else if (init == "random") {
      embedding <- rand_init(n_vertices, n_components)
    }
    else if (init == "normlaplacian") {
      embedding <- normalized_laplacian_init(V, ndim = n_components,
                                             use_RSpectra = TRUE,
                                             verbose = verbose)
    }
    else if (init == "laplacian") {
      embedding <- laplacian_eigenmap(V, ndim = n_components, use_RSpectra = TRUE,
                                      verbose = verbose)
    }
    else {
      embedding <- scaled_pca(X, ndim = n_components, verbose = verbose)
    }
  }

  if (is.null(n_epochs) || n_epochs <= 0) {
    if (n_vertices <= 10000) {
      n_epochs <- 500
    }
    else {
      n_epochs <- 200
    }
  }

  V@x[V@x < max(V@x) / n_epochs] <- 0
  V <- Matrix::drop0(V)
  epochs_per_sample <- make_epochs_per_sample(V@x, n_epochs)

  positive_head <- V@i
  positive_tail <- Matrix::which(V != 0, arr.ind = TRUE)[, 2] - 1


  optimize_layout_cpp(embedding, positive_head, positive_tail,
                  n_epochs, n_vertices,
                  epochs_per_sample, a, b, gamma,
                  initial_alpha = alpha, negative_sample_rate,
                  verbose = verbose)
  embedding
}


make_epochs_per_sample <- function(weights, n_epochs) {
  result <- rep(-1, length(weights))
  n_samples = n_epochs * (weights / max(weights))
  result[n_samples > 0] <- n_epochs / n_samples[n_samples > 0]
  result
}

find_ab_params <- function(spread = 1, min_dist = 0.001) {
  xv <- seq(from = 0, to = spread * 3, length.out = 300)
  yv <- rep(0, length(xv))
  yv[xv < min_dist] <- 1
  yv[xv >= min_dist] <- exp(-(xv[xv >= min_dist] - min_dist) / spread)
  result <- try({
    stats::nls(yv ~ 1 / (1 + a * xv ^ (2 * b)),
               start = list(a = 1, b = 1))$m$getPars()
  }, silent = TRUE)
  if (class(result) == "try-error") {
    stop("Can't find a, b for provided spread/min_dist values")
  }
  result
}

