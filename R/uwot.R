#' @useDynLib uwot
#' @importFrom Rcpp sourceCpp
NULL

#' Dimensionality Reduction with UMAP
#'
#' Carry out dimensionality reduction of a dataset using the Uniform Manifold
#' Approximation and Projection (UMAP) method (McInnes & Healy, 2018). Some of
#' the following help text is lifted verbatim from the Python reference
#' implementation at \url{https://github.com/lmcinnes/umap}.
#'
#' @param X Input data. Can be a \code{data.frame}, \code{matrix} or \code{dist}
#'   object. Matrix and data frames should contain one observation per row. Data
#'   frames will have any non-numeric columns removed.
#' @param n_neighbors The size of local neighborhood (in terms of number of
#'   neighboring sample points) used for manifold approximation. Larger values
#'   result in more global views of the manifold, while smaller values result in
#'   more local data being preserved. In general values should be in the range
#'   \code{2} to \code{100}.
#' @param n_components The dimension of the space to embed into. This defaults
#'   to \code{2} to provide easy visualization, but can reasonably be set to any
#'   integer value in the range \code{2} to \code{100}.
#' @param n_epochs Number of epochs to use during the optimization of the
#'   embedded coordinates. By default, this value is set to \code{500} for datasets
#'   containing 10,000 vertices or less, and \code{200} otherwise.
#' @param alpha Initial learning rate used in optimization of the coordinates.
#' @param init Type of initialization for the coordinates. Options are:
#'   \itemize{
#'     \item \code{"spectral"} Spectral embedding using the normalized Laplacian
#'     of the fuzzy 1-skeleton, with Gaussian noise added.
#'     \item \code{"normlaplacian"}. Spectral embedding using the normalized
#'     Laplacian of the fuzzy 1-skeleton, without noise.
#'     \item \code{random}. Coordinates assigned using a uniform random
#'     distribution between -10 and 10.
#'     \item \code{"laplacian"}. Spectral embedding using the Laplacian Eigenmap
#'     (Belkin and Niyogi, 2002).
#'     \item \code{"spca"}. The first two principal components from PCA of
#'     \code{X} if \code{X} is a data frame, and from a 2-dimensional classical
#'     MDS if \code{X} is of class \code{"dist"}. These vectors are then scaled
#'     so their standard deviation is 0.0001.
#'     \item A matrix of initial coordinates.
#'   }
#' @param spread The effective scale of embedded points. In combination with
#'   \code{min_dist}, this determines how clustered/clumped the embedded points
#'   are.
#' @param min_dist The effective minimum distance between embedded points.
#'   Smaller values will result in a more clustered/clumped embedding where
#'   nearby points on the manifold are drawn closer together, while larger
#'   values will result on a more even dispersal of points. The value should be
#'   set relative to the \code{spread} value, which determines the scale at
#'   which embedded points will be spread out.
#' @param set_op_mix_ratio Interpolate between (fuzzy) union and intersection as
#'   the set operation used to combine local fuzzy simplicial sets to obtain a
#'   global fuzzy simplicial sets. Both fuzzy set operations use the product
#'   t-norm. The value of this parameter should be between \code{0.0} and
#'   \code{1.0}; a value of \code{1.0} will use a pure fuzzy union, while
#'   \code{0.0} will use a pure fuzzy intersection.
#' @param local_connectivity The local connectivity required -- i.e. the number
#'   of nearest neighbors that should be assumed to be connected at a local
#'   level. The higher this value the more connected the manifold becomes
#'   locally. In practice this should be not more than the local intrinsic
#'   dimension of the manifold.
#' @param bandwidth The effective bandwidth of the kernel if we view the
#'   algorithm as similar to Laplacian eigenmaps. Larger values induce more
#'   connectivity and a more global view of the data, smaller values concentrate
#'   more locally.
#' @param gamma Weighting applied to negative samples in low dimensional embedding
#'   optimization. Values higher than one will result in greater weight
#'   being given to negative samples.
#' @param negative_sample_rate The number of negative edge/1-simplex samples to
#'   use per positive edge/1-simplex sample in optimizing the low dimensional
#'   embedding.
#' @param a More specific parameters controlling the embedding. If \code{NULL}
#'   these values are set automatically as determined by \code{min_dist} and
#'   \code{spread}.
#' @param b More specific parameters controlling the embedding. If \code{NULL}
#'   these values are set automatically as determined by \code{min_dist} and
#'   \code{spread}.
#' @param nn_method Method for finding nearest neighbors. Options are:
#'   \itemize{
#'     \item \code{"fnn"}. Use exact nearest neighbors via the
#'       \href{https://cran.r-project.org/package=FNN}{FNN} package.
#'     \item \code{"annoy"} Use approximate nearest neighbors via the
#'       \href{https://cran.r-project.org/package=RcppAnnoy}{RcppAnnoy} package.
#'    }
#'   By default, if \code{X} has less than 4,096 vertices, the exact nearest
#'   neighbors are found. Otherwise, approximate nearest neighbors are used.
#' @param verbose If \code{TRUE}, log details to the console.
#' @return A matrix of optimized coordinates.
#' @examples
#' \dontrun{
#' iris_umap <- umap(iris, n_neighbors = 50, alpha = 0.5, init = "random")
#'
#' # Load mnist from somewhere, e.g.
#' # devtools::install_github("jlmelville/snedata")
#' # mnist <- snedata::download_mnist()
#' mnist_umap <- umap(mnist, n_neighbors = 15, min_dist = 0.001, verbose = TRUE)
#' }
#' @references
#' Belkin, M., & Niyogi, P. (2002).
#' Laplacian eigenmaps and spectral techniques for embedding and clustering.
#' In \emph{Advances in neural information processing systems}
#' (pp. 585-591).
#' \url{http://papers.nips.cc/paper/1961-laplacian-eigenmaps-and-spectral-techniques-for-embedding-and-clustering.pdf}
#'
#' McInnes, L., & Healey, J. (2018).
#' UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction
#' \emph{arXiv preprint} \emph{arXiv}:1802.03426.
#' \url{https://arxiv.org/abs/1802.03426}
#'
#' @export
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

  if (n_neighbors < 2) {
    stop("n_neighbors must be >= 2")
  }

  if (set_op_mix_ratio < 0.0 || set_op_mix_ratio > 1.0) {
    stop("set_op_mix_ratio must be between 0.0 and 1.0")
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
      embedding <- spectral_init(V, ndim = n_components, verbose = verbose)
    }
    else if (init == "random") {
      embedding <- rand_init(n_vertices, n_components)
    }
    else if (init == "normlaplacian") {
      embedding <- normalized_laplacian_init(V, ndim = n_components,
                                             verbose = verbose)
    }
    else if (init == "laplacian") {
      embedding <- laplacian_eigenmap(V, ndim = n_components, verbose = verbose)
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

