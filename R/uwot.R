#' Dimensionality Reduction with UMAP
#'
#' Carry out dimensionality reduction of a dataset using the Uniform Manifold
#' Approximation and Projection (UMAP) method (McInnes & Healy, 2018). Some of
#' the following help text is lifted verbatim from the Python reference
#' implementation at \url{https://github.com/lmcinnes/umap}.
#'
#' @param X Input data. Can be a \code{\link{data.frame}}, \code{\link{matrix}},
#'   \code{\link[stats]{dist}} object or \code{\link[Matrix]{sparseMatrix}}.
#'   A sparse matrix is interpreted as a distance matrix and both implicit and
#'   explicit zero entries are ignored. Set zero distances you want to keep to
#'   an arbitrarily small non-zero value (e.g. \code{1e-10}). Matrix and data
#'   frames should contain one observation per row. Data frames will have any
#'   non-numeric columns removed.
#' @param n_neighbors The size of local neighborhood (in terms of number of
#'   neighboring sample points) used for manifold approximation. Larger values
#'   result in more global views of the manifold, while smaller values result in
#'   more local data being preserved. In general values should be in the range
#'   \code{2} to \code{100}.
#' @param n_components The dimension of the space to embed into. This defaults
#'   to \code{2} to provide easy visualization, but can reasonably be set to any
#'   integer value in the range \code{2} to \code{100}.
#' @param metric Type of distance metric to use to find nearest neighbors. One
#'   of:
#' \itemize{
#'   \item \code{"euclidean"} (the default)
#'   \item \code{"cosine"}
#'   \item \code{"manhattan"}
#' }
#' Only applies if \code{nn_method = "annoy"} (for \code{nn_method = "fnn"}, the
#' distance metric is always "euclidean").
#' @param n_epochs Number of epochs to use during the optimization of the
#'   embedded coordinates. By default, this value is set to \code{500} for datasets
#'   containing 10,000 vertices or less, and \code{200} otherwise.
#' @param scale Scaling to apply to \code{X} if it is a data frame or matrix:
#' \itemize{
#'   \item{\code{"none"} or \code{FALSE} or \code{NULL}} No scaling.
#'   \item{\code{"scale"} or \code{TRUE}} Scale each column to zero mean and variance 1.
#'   \item{\code{"maxabs"}} Center each column to mean 0, then divide each element by the
#'   maximum absolute value over the entire matrix.
#'   \item{\code{"range"}} Range scale the entire matrix, so the smallest element is 0 and
#'   the largest is 1.
#' }
#' For UMAP, the default is \code{"none"}.
#' @param alpha Initial learning rate used in optimization of the coordinates.
#' @param init Type of initialization for the coordinates. Options are:
#'   \itemize{
#'     \item \code{"spectral"} Spectral embedding using the normalized Laplacian
#'     of the fuzzy 1-skeleton, with Gaussian noise added.
#'     \item \code{"normlaplacian"}. Spectral embedding using the normalized
#'     Laplacian of the fuzzy 1-skeleton, without noise.
#'     \item \code{"random"}. Coordinates assigned using a uniform random
#'     distribution between -10 and 10.
#'     \item \code{"lvrandom"}. Coordinates assigned using a Gaussian
#'     distribution with standard deviation 1e-4, as used in LargeVis
#'     (Tang et al., 2016) and t-SNE.
#'     \item \code{"laplacian"}. Spectral embedding using the Laplacian Eigenmap
#'     (Belkin and Niyogi, 2002).
#'     \item \code{"pca"}. The first two principal components from PCA of
#'     \code{X} if \code{X} is a data frame, and from a 2-dimensional classical
#'     MDS if \code{X} is of class \code{"dist"}.
#'     \item \code{"spca"}. Like \code{"pca"}, but each dimension is then scaled
#'     so the standard deviation is 1e-4, to give a distribution similar to
#'     that used in t-SNE.
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
#' @param n_trees Number of trees to build when constructing the nearest
#'   neighbor index. The more trees specified, the larger the index, but the
#'   better the results. With \code{search_k}, determines the accuracy of the
#'   Annoy nearest neighbor search. Only used if the \code{nn_method} is
#'   \code{"annoy"}. Sensible values are between \code{10} to \code{100}.
#' @param search_k Number of nodes to search during the neighbor retrieval. The
#'   larger k, the more the accurate results, but the longer the search takes.
#'   With \code{n_trees}, determines the accuracy of the Annoy nearest neighbor
#'   search. Only used if the \code{nn_method} is \code{"annoy"}.
#' @param approx_pow If \code{TRUE}, use an approximation to the power function
#'   in the UMAP gradient, from
#'   \url{https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/}.
#' @param y Optional target array for supervised dimension reduction. Must be a
#'   factor or numeric vector with the same length as \code{X}.
#' @param target_n_neighbors Number of nearest neighbors to use to construct the
#'   target simplcial set. Default value is \code{n_neighbors}. Applies only if
#'   \code{y} is non-\code{NULL} and \code{numeric}.
#' @param target_weight Weighting factor between data topology and target
#'   topology. A value of 0.0 weights entirely on data, a value of 1.0 weights
#'   entirely on target. The default of 0.5 balances the weighting equally
#'   between data and target. Only applies if \code{y} is non-\code{NULL}.
#' @param ret_model If \code{TRUE}, then return extra data that can be used to
#'   add new data to an existing embedding via \code{\link{umap_transform}}.
#'   Otherwise, just return the coordinates.
#' @param n_threads Number of threads to use. Default is half that recommended
#'   by RcppParallel. For nearest neighbor search, only applies if
#'   \code{nn_method = "annoy"}.
#' @param grain_size Minimum batch size for multithreading. If the number of
#'   items to process in a thread falls below this number, then no threads will
#'   be used. Used in conjunction with \code{n_threads}.
#' @param verbose If \code{TRUE}, log details to the console.
#' @return A matrix of optimized coordinates, or if \code{ret_model = TRUE}, a
#'   list containing extra information that can be used to add new data to an
#'   existing embedding via \code{\link{umap_transform}}. In this case, the
#'   coordinates are available in the list item \code{embedding}.
#' @examples
#' \dontrun{
#' iris_umap <- umap(iris, n_neighbors = 50, alpha = 0.5, init = "random")
#'
#' # Faster approximation to the gradient
#' iris_umap <- umap(iris, n_neighbors = 15, approx_pow = TRUE)
#'
#' # Load mnist from somewhere, e.g.
#' # devtools::install_github("jlmelville/snedata")
#' # mnist <- snedata::download_mnist()
#' mnist_umap <- umap(mnist, n_neighbors = 15, min_dist = 0.001, verbose = TRUE)
#'
#' # Supervised dimension reduction
#' mnist_sumap <- umap(mnist, n_neighbors = 15, min_dist = 0.001, verbose = TRUE,
#'                     y = mnist$Label, target_weight = 0.5)
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
#' Tang, J., Liu, J., Zhang, M., & Mei, Q. (2016, April).
#' Visualizing large-scale and high-dimensional data.
#' In \emph{Proceedings of the 25th International Conference on World Wide Web}
#' (pp. 287-297).
#' International World Wide Web Conferences Steering Committee.
#' \url{https://arxiv.org/abs/1602.00370}
#'
#' Van der Maaten, L., & Hinton, G. (2008).
#' Visualizing data using t-SNE.
#' \emph{Journal of Machine Learning Research}, \emph{9} (2579-2605).
#' \url{http://www.jmlr.org/papers/v9/vandermaaten08a.html}
#' @export
umap <- function(X, n_neighbors = 15, n_components = 2, metric = "euclidean",
                 n_epochs = NULL, alpha = 1, scale = FALSE, init = "spectral",
                 spread = 1, min_dist = 0.01,
                 set_op_mix_ratio = 1.0, local_connectivity = 1.0,
                 bandwidth = 1.0, gamma = 1.0,
                 negative_sample_rate = 5.0, a = NULL, b = NULL,
                 nn_method = NULL, n_trees = 50,
                 search_k = 2 * n_neighbors * n_trees,
                 approx_pow = FALSE,
                 y = NULL, target_n_neighbors = n_neighbors,
                 target_weight = 0.5,
                 ret_model = FALSE,
                 n_threads = max(1, RcppParallel::defaultNumThreads() / 2),
                 grain_size = 1,
                 verbose = getOption("verbose", TRUE)) {
  uwot(
    X = X, n_neighbors = n_neighbors, n_components = n_components,
    metric = metric, n_epochs = n_epochs, alpha = alpha, scale = scale,
    init = init, spread = spread, min_dist = min_dist,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity, bandwidth = bandwidth,
    gamma = gamma, negative_sample_rate = negative_sample_rate,
    a = a, b = b, nn_method = nn_method, n_trees = n_trees,
    search_k = search_k, method = "umap", approx_pow = approx_pow,
    n_threads = n_threads, grain_size = grain_size,
    y = y, target_n_neighbors = target_n_neighbors,
    target_weight = target_weight,
    ret_model = ret_model,
    verbose = verbose
  )
}

#' Dimensionality Reduction Using t-Distributed UMAP (t-UMAP)
#'
#' A faster (but less flexible) version of the UMAP gradient. FOr more detail on
#' UMAP, see the  \code{\link{umap}} function.
#'
#' By setting the UMAP curve parameters \code{a} and \code{b} to \code{1}, you
#' get back the Cauchy distribution as used in t-SNE and LargeVis. It also
#' results in a substantially simplified gradient expression. This can give
#' a speed improvement of around 50\%.
#'
#' @param X Input data. Can be a \code{\link{data.frame}}, \code{\link{matrix}},
#'   \code{\link[stats]{dist}} object or \code{\link[Matrix]{sparseMatrix}}.
#'   A sparse matrix is interpreted as a distance matrix and both implicit and
#'   explicit zero entries are ignored. Set zero distances you want to keep to
#'   an arbitrarily small non-zero value (e.g. \code{1e-10}). Matrix and data
#'   frames should contain one observation per row. Data frames will have any
#'   non-numeric columns removed.
#' @param n_neighbors The size of local neighborhood (in terms of number of
#'   neighboring sample points) used for manifold approximation. Larger values
#'   result in more global views of the manifold, while smaller values result in
#'   more local data being preserved. In general values should be in the range
#'   \code{2} to \code{100}.
#' @param n_components The dimension of the space to embed into. This defaults
#'   to \code{2} to provide easy visualization, but can reasonably be set to any
#'   integer value in the range \code{2} to \code{100}.
#' @param metric Type of distance metric to use to find nearest neighbors. One
#'   of:
#' \itemize{
#'   \item \code{"euclidean"} (the default)
#'   \item \code{"cosine"}
#'   \item \code{"manhattan"}
#' }
#' Only applies if \code{nn_method = "annoy"} (for \code{nn_method = "fnn"}, the
#' distance metric is always "euclidean").
#' @param n_epochs Number of epochs to use during the optimization of the
#'   embedded coordinates. By default, this value is set to \code{500} for datasets
#'   containing 10,000 vertices or less, and \code{200} otherwise.
#' @param alpha Initial learning rate used in optimization of the coordinates.
#' @param scale Scaling to apply to \code{X} if it is a data frame or matrix:
#' \itemize{
#'   \item{\code{"none"} or \code{FALSE} or \code{NULL}} No scaling.
#'   \item{\code{"scale"} or \code{TRUE}} Scale each column to zero mean and variance 1.
#'   \item{\code{"maxabs"}} Center each column to mean 0, then divide each element by the
#'   maximum absolute value over the entire matrix.
#'   \item{\code{"range"}} Range scale the entire matrix, so the smallest element is 0 and
#'   the largest is 1.
#' }
#' For t-UMAP, the default is \code{"none"}.
#' @param init Type of initialization for the coordinates. Options are:
#'   \itemize{
#'     \item \code{"spectral"} Spectral embedding using the normalized Laplacian
#'     of the fuzzy 1-skeleton, with Gaussian noise added.
#'     \item \code{"normlaplacian"}. Spectral embedding using the normalized
#'     Laplacian of the fuzzy 1-skeleton, without noise.
#'     \item \code{"random"}. Coordinates assigned using a uniform random
#'     distribution between -10 and 10.
#'     \item \code{"lvrandom"}. Coordinates assigned using a Gaussian
#'     distribution with standard deviation 1e-4, as used in LargeVis
#'     (Tang et al., 2016) and t-SNE.
#'     \item \code{"laplacian"}. Spectral embedding using the Laplacian Eigenmap
#'     (Belkin and Niyogi, 2002).
#'     \item \code{"pca"}. The first two principal components from PCA of
#'     \code{X} if \code{X} is a data frame, and from a 2-dimensional classical
#'     MDS if \code{X} is of class \code{"dist"}.
#'     \item \code{"spca"}. Like \code{"pca"}, but each dimension is then scaled
#'     so the standard deviation is 1e-4, to give a distribution similar to
#'     that used in t-SNE.
#'     \item A matrix of initial coordinates.
#'   }
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
#' @param nn_method Method for finding nearest neighbors. Options are:
#'   \itemize{
#'     \item \code{"fnn"}. Use exact nearest neighbors via the
#'       \href{https://cran.r-project.org/package=FNN}{FNN} package.
#'     \item \code{"annoy"} Use approximate nearest neighbors via the
#'       \href{https://cran.r-project.org/package=RcppAnnoy}{RcppAnnoy} package.
#'    }
#'   By default, if \code{X} has less than 4,096 vertices, the exact nearest
#'   neighbors are found. Otherwise, approximate nearest neighbors are used.
#' @param n_trees Number of trees to build when constructing the nearest
#'   neighbor index. The more trees specified, the larger the index, but the
#'   better the results. With \code{search_k}, determines the accuracy of the
#'   Annoy nearest neighbor search. Only used if the \code{nn_method} is
#'   \code{"annoy"}. Sensible values are between \code{10} to \code{100}.
#' @param search_k Number of nodes to search during the neighbor retrieval. The
#'   larger k, the more the accurate results, but the longer the search takes.
#'   With \code{n_trees}, determines the accuracy of the Annoy nearest neighbor
#'   search. Only used if the \code{nn_method} is \code{"annoy"}.
#' @param y Optional target array for supervised dimension reduction. Must be a
#'   factor or numeric vector with the same length as \code{X}.
#' @param target_n_neighbors Number of nearest neighbors to use to construct the
#'   target simplcial set. Default value is \code{n_neighbors}. Applies only if
#'   \code{y} is non-\code{NULL} and \code{numeric}.
#' @param target_weight Weighting factor between data topology and target
#'   topology. A value of 0.0 weights entirely on data, a value of 1.0 weights
#'   entirely on target. The default of 0.5 balances the weighting equally
#'   between data and target. Only applies if \code{y} is non-\code{NULL}.
#' @param ret_model If \code{TRUE}, then return extra data that can be used to
#'   add new data to an existing embedding via \code{\link{umap_transform}}.
#'   Otherwise, just return the coordinates.
#' @param n_threads Number of threads to use. Default is half that recommended
#'   by RcppParallel. For nearest neighbor search, only applies if
#'   \code{nn_method = "annoy"}.
#' @param grain_size Minimum batch size for multithreading. If the number of
#'   items to process in a thread falls below this number, then no threads will
#'   be used. Used in conjunction with \code{n_threads}.
#' @param verbose If \code{TRUE}, log details to the console.
#' @return A matrix of optimized coordinates, or if \code{ret_model = TRUE}, a
#'   list containing extra information that can be used to add new data to an
#'   existing embedding via \code{\link{umap_transform}}. In this case, the
#'   coordinates are available in the list item \code{embedding}.
#' @export
tumap <- function(X, n_neighbors = 15, n_components = 2, metric = "euclidean",
                  n_epochs = NULL,
                  alpha = 1, scale = FALSE, init = "spectral",
                  set_op_mix_ratio = 1.0, local_connectivity = 1.0,
                  bandwidth = 1.0, gamma = 1.0,
                  negative_sample_rate = 5.0,
                  nn_method = NULL, n_trees = 50,
                  search_k = 2 * n_neighbors * n_trees,
                  n_threads = max(1, RcppParallel::defaultNumThreads() / 2),
                  grain_size = 1,
                  y = NULL, target_n_neighbors = n_neighbors,
                  target_weight = 0.5,
                  ret_model = FALSE,
                  verbose = getOption("verbose", TRUE)) {
  uwot(
    X = X, n_neighbors = n_neighbors, n_components = n_components,
    metric = metric,
    n_epochs = n_epochs, alpha = alpha, scale = scale, init = init,
    spread = NULL, min_dist = NULL, set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity, bandwidth = bandwidth,
    gamma = gamma, negative_sample_rate = negative_sample_rate,
    a = NULL, b = NULL, nn_method = nn_method, n_trees = n_trees,
    search_k = search_k, method = "tumap",
    n_threads = n_threads, grain_size = grain_size,
    y = y, target_n_neighbors = target_n_neighbors,
    target_weight = target_weight,
    ret_model = ret_model,
    verbose = verbose
  )
}

#' Dimensionality Reduction with a LargeVis-like method
#'
#' Carry out dimensionality reduction of a dataset using a method similar to
#' LargeVis (Tang et al., 2016).
#'
#' \code{lvish} differs from the official LargeVis implementation in the
#' following:
#'
#' \itemize{
#'   \item Only the nearest-neighbor index search phase is multi-threaded.
#'   \item Matrix input data is not normalized.
#'   \item The \code{n_trees} parameter cannot be dynamically chosen based on
#'   data set size.
#'   \item Nearest neighbor results are not refined via the
#'   neighbor-of-my-neighbor method. The \code{search_k} parameter is twice
#'   as large than default to compensate.
#'   \item Gradient values are clipped to \code{4.0} rather than \code{5.0}.
#'   \item Negative edges are generated by uniform sampling of vertexes rather
#'   than their degree ^ 0.75.
#'   \item The default number of samples is much reduced. The default number of
#'   epochs, \code{n_epochs}, is set to \code{5000}, much larger than for
#'   \code{\link{umap}}, but may need to be increased further depending on your
#'   dataset. Using \code{init = "spectral"} can help.
#' }
#'
#' @param X Input data. Can be a \code{\link{data.frame}}, \code{\link{matrix}},
#'   \code{\link[stats]{dist}} object or \code{\link[Matrix]{sparseMatrix}}.
#'   A sparse matrix is interpreted as a distance matrix and both implicit and
#'   explicit zero entries are ignored. Set zero distances you want to keep to
#'   an arbitrarily small non-zero value (e.g. \code{1e-10}). Matrix and data
#'   frames should contain one observation per row. Data frames will have any
#'   non-numeric columns removed.
#' @param perplexity Controls the size of the local neighborhood used for
#'   manifold approximation. This is the analogous to \code{n_neighbors} in
#'   \code{\link{umap}}. Change this, rather than \code{n_neighbors}.
#' @param n_neighbors The number of neighbors to use when calculating the
#'  \code{perplexity}. Usually set to three times the value of the
#'  \code{perplexity}. Must be at least as large as \code{perplexity}.
#' @param n_components The dimension of the space to embed into. This defaults
#'   to \code{2} to provide easy visualization, but can reasonably be set to any
#'   integer value in the range \code{2} to \code{100}.
#' @param metric Type of distance metric to use to find nearest neighbors. One
#'   of:
#' \itemize{
#'   \item \code{"euclidean"} (the default)
#'   \item \code{"cosine"}
#'   \item \code{"manhattan"}
#' }
#' Only applies if \code{nn_method = "annoy"} (for \code{nn_method = "fnn"}, the
#' distance metric is always "euclidean").
#' @param n_epochs Number of epochs to use during the optimization of the
#'   embedded coordinates. The default is calculate the number of epochs
#'   dynamically based on dataset size, to give the same number of edge samples
#'   as the LargeVis defaults. This is usually substantially larger than the
#'   UMAP defaults.
#' @param alpha Initial learning rate used in optimization of the coordinates.
#' @param scale Scaling to apply to \code{X} if it is a data frame or matrix:
#' \itemize{
#'   \item{\code{"none"} or \code{FALSE} or \code{NULL}} No scaling.
#'   \item{\code{"scale"} or \code{TRUE}} Scale each column to zero mean and variance 1.
#'   \item{\code{"maxabs"}} Center each column to mean 0, then divide each element by the
#'   maximum absolute value over the entire matrix.
#'   \item{\code{"range"}} Range scale the entire matrix, so the smallest element is 0 and
#'   the largest is 1.
#' }
#' For lvish, the default is \code{"maxabs"}, for consistency with LargeVis.
#' @param init Type of initialization for the coordinates. Options are:
#'   \itemize{
#'     \item \code{"spectral"} Spectral embedding using the normalized Laplacian
#'     of the fuzzy 1-skeleton, with Gaussian noise added.
#'     \item \code{"normlaplacian"}. Spectral embedding using the normalized
#'     Laplacian of the fuzzy 1-skeleton, without noise.
#'     \item \code{"random"}. Coordinates assigned using a uniform random
#'     distribution between -10 and 10.
#'     \item \code{"lvrandom"}. Coordinates assigned using a Gaussian
#'     distribution with standard deviation 1e-4, as used in LargeVis
#'     (Tang et al., 2016) and t-SNE.
#'     \item \code{"laplacian"}. Spectral embedding using the Laplacian Eigenmap
#'     (Belkin and Niyogi, 2002).
#'     \item \code{"pca"}. The first two principal components from PCA of
#'     \code{X} if \code{X} is a data frame, and from a 2-dimensional classical
#'     MDS if \code{X} is of class \code{"dist"}.
#'     \item \code{"spca"}. Like \code{"pca"}, but each dimension is then scaled
#'     so the standard deviation is 1e-4, to give a distribution similar to
#'     that used in t-SNE and LargeVis.
#'     \item A matrix of initial coordinates.
#'   }
#' @param gamma Weighting applied to negative samples in low dimensional embedding
#'   optimization. Values higher than one will result in greater weight
#'   being given to negative samples.
#' @param negative_sample_rate The number of negative edge/1-simplex samples to
#'   use per positive edge/1-simplex sample in optimizing the low dimensional
#'   embedding.
#' @param nn_method Method for finding nearest neighbors. Options are:
#'   \itemize{
#'     \item \code{"fnn"}. Use exact nearest neighbors via the
#'       \href{https://cran.r-project.org/package=FNN}{FNN} package.
#'     \item \code{"annoy"} Use approximate nearest neighbors via the
#'       \href{https://cran.r-project.org/package=RcppAnnoy}{RcppAnnoy} package.
#'    }
#'   By default, if \code{X} has less than 4,096 vertices, the exact nearest
#'   neighbors are found. Otherwise, approximate nearest neighbors are used.
#' @param n_trees Number of trees to build when constructing the nearest
#'   neighbor index. The more trees specified, the larger the index, but the
#'   better the results. With \code{search_k}, determines the accuracy of the
#'   Annoy nearest neighbor search. Only used if the \code{nn_method} is
#'   \code{"annoy"}. Sensible values are between \code{10} to \code{100}.
#' @param search_k Number of nodes to search during the neighbor retrieval. The
#'   larger k, the more the accurate results, but the longer the search takes.
#'   With \code{n_trees}, determines the accuracy of the Annoy nearest neighbor
#'   search. Only used if the \code{nn_method} is \code{"annoy"}.
#' @param n_threads Number of threads to use. Default is half that recommended
#'   by RcppParallel. For nearest neighbor search, only applies if
#'   \code{nn_method = "annoy"}.
#' @param grain_size Minimum batch size for multithreading. If the number of
#'   items to process in a thread falls below this number, then no threads will
#'   be used. Used in conjunction with \code{n_threads}.
#' @param kernel Type of kernel function to create input probabilities. Can be
#'   one of \code{"gauss"} (the default) or \code{"knn"}. \code{"gauss"} uses
#'   the usual Gaussian weighted similarities. \code{"knn"} assigns equal
#'   probabilities to every edge in the nearest neighbor graph, and zero
#'   otherwise, using \code{perplexity} nearest neighbors. The \code{n_neighbors}
#'   parameter is ignored in this case.
#' @param verbose If \code{TRUE}, log details to the console.
#' @return A matrix of optimized coordinates.
#' @references
#' Tang, J., Liu, J., Zhang, M., & Mei, Q. (2016, April).
#' Visualizing large-scale and high-dimensional data.
#' In \emph{Proceedings of the 25th International Conference on World Wide Web}
#' (pp. 287-297).
#' International World Wide Web Conferences Steering Committee.
#' \url{https://arxiv.org/abs/1602.00370}
#' @examples
#' \dontrun{
#' # Use perplexity rather than n_neighbors to control the size of the local
#' neighborhood iris_lvish <- umap(iris, perplexity = 50, alpha = 0.5,
#'                                 init = "random")
#'
#' # Default number of epochs is much larger than for UMAP, assumes random
#' # initialization
#' # If using a more global initialization, can use fewer epochs
#' iris_lvish_short <- umap(iris, perpelxity = 50, n_epochs = 1000)
#' }
#' @export
lvish <- function(X, perplexity = 50, n_neighbors = perplexity * 3,
                  n_components = 2, metric = "euclidean", n_epochs = -1,
                  alpha = 1, scale = "maxabs", init = "lvrandom", gamma = 7,
                  negative_sample_rate = 5.0,
                  nn_method = NULL, n_trees = 50,
                  search_k = 2 * n_neighbors * n_trees,
                  n_threads = max(1, RcppParallel::defaultNumThreads() / 2),
                  grain_size = 1,
                  kernel = "gauss",
                  verbose = getOption("verbose", TRUE)) {
  uwot(X,
    n_neighbors = n_neighbors, n_components = n_components,
    metric = metric,
    n_epochs = n_epochs, alpha = alpha, scale = scale, init = init,
    gamma = gamma, negative_sample_rate = negative_sample_rate,
    nn_method = nn_method, n_trees = n_trees, search_k = search_k,
    method = "largevis", perplexity = perplexity,
    n_threads = n_threads,
    grain_size = grain_size, kernel = kernel, verbose = verbose
  )
}

# Function that does all the real work
uwot <- function(X, n_neighbors = 15, n_components = 2, metric = "euclidean",
                 n_epochs = NULL,
                 alpha = 1, scale = FALSE, init = "spectral",
                 spread = 1, min_dist = 0.01,
                 set_op_mix_ratio = 1.0, local_connectivity = 1.0,
                 bandwidth = 1.0, gamma = 1.0,
                 negative_sample_rate = 5.0, a = NULL, b = NULL,
                 nn_method = NULL, n_trees = 50,
                 search_k = 2 * n_neighbors * n_trees,
                 method = "umap", perplexity = 50, approx_pow = FALSE,
                 y = NULL, target_n_neighbors = n_neighbors,
                 target_weight = 0.5,
                 n_threads = max(1, RcppParallel::defaultNumThreads() / 2),
                 kernel = "gauss",
                 grain_size = 1,
                 ret_model = FALSE,
                 verbose = getOption("verbose", TRUE)) {
  if (method == "umap" && (is.null(a) || is.null(b))) {
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

  if (local_connectivity < 1.0) {
    stop("local_connectivity cannot be < 1.0")
  }

  if (n_threads > 0) {
    RcppParallel::setThreadOptions(numThreads = n_threads)
  }

  if (methods::is(X, "dist")) {
    if (ret_model) {
      stop("Can only create models with dense matrix or data frame input")
    }
    n_vertices <- attr(X, "Size")
    tsmessage("Read ", n_vertices, " rows")
  }
  else if (methods::is(X, "sparseMatrix")) {
    if (ret_model) {
      stop("Can only create models with dense matrix or data frame input")
    }
    n_vertices <- nrow(X)
    if (ncol(X) != n_vertices) {
      stop("Sparse matrices are only supported as distance matrices")
    }
    tsmessage("Read ", n_vertices, " rows of sparse distance matrix")
  }
  else {
    if (methods::is(X, "data.frame")) {
      indexes <- which(vapply(X, is.numeric, logical(1)))
      if (length(indexes) == 0) {
        stop("No numeric columns found")
      }
      X <- as.matrix(X[, indexes])
    }
    n_vertices <- nrow(X)
    tsmessage(
      "Read ", n_vertices, " rows and found ", ncol(X),
      " numeric columns"
    )
    X <- scale_input(X,
      scale_type = scale, ret_model = ret_model,
      verbose = verbose
    )
  }

  if (method == "largevis" && kernel == "knn") {
    n_neighbors <- perplexity
  }

  if (n_neighbors > n_vertices) {
    # for LargeVis, n_neighbors normally determined from perplexity
    # not an error to be too large
    if (method == "largevis") {
      tsmessage("Setting n_neighbors to ", n_vertices)
      n_neighbors <- n_vertices
    }
    else {
      stop("n_neighbors must be smaller than the dataset size")
    }
  }

  metric <- match.arg(tolower(metric), c("euclidean", "cosine", "manhattan"))

  if (is.null(nn_method)) {
    if (n_vertices < 4096 && metric == "euclidean" && !ret_model) {
      tsmessage("Using FNN for neighbor search, n_neighbors = ", n_neighbors)
      nn_method <- "fnn"
    }
    else {
      tsmessage("Using Annoy for neighbor search, n_neighbors = ", n_neighbors)
      nn_method <- "annoy"
    }
  }
  nn_method <- match.arg(tolower(nn_method), c("annoy", "fnn"))
  if (nn_method == "fnn" && metric != "euclidean") {
    stop(
      "nn_method = 'FNN' is only compatible with distance metric ",
      "'euclidean'"
    )
  }
  if (nn_method == "fnn" && ret_model) {
    stop("nn_method = 'FNN' is incompatible with ret_model = TRUE")
  }
  nn <- find_nn(X, n_neighbors,
    method = nn_method, metric = metric,
    n_trees = n_trees,
    n_threads = n_threads, grain_size = grain_size,
    search_k = search_k, ret_index = ret_model, verbose = verbose
  )
  gc()
  if (any(is.infinite(nn$dist))) {
    stop("Infinite distances found in nearest neighbors")
  }

  if (method == "largevis") {
    if (perplexity >= n_vertices) {
      stop("perplexity can be no larger than ", n_vertices - 1)
    }
    V <- perplexity_similarities(
      nn = nn, perplexity = perplexity,
      n_threads = n_threads,
      grain_size = grain_size,
      kernel = kernel,
      verbose = verbose
    )
  }
  else {
    V <- fuzzy_simplicial_set(
      nn = nn,
      set_op_mix_ratio = set_op_mix_ratio,
      local_connectivity = local_connectivity,
      bandwidth = bandwidth,
      n_threads = n_threads,
      grain_size = grain_size,
      verbose = verbose
    )
  }
  if (any(is.na(V))) {
    stop("Non-finite entries in the input matrix")
  }
  gc()

  if (!is.null(y)) {
    if (is.factor(y)) {
      if (target_weight < 1.0) {
        far_dist <- 2.5 * (1.0 / (1.0 - target_weight))
      }
      else {
        far_dist <- 1.0e12
      }
      tsmessage(
        "Applying categorical set intersection, target weight = ",
        formatC(target_weight), " far distance = ", formatC(far_dist)
      )

      V <- categorical_simplicial_set_intersection(V, y,
        far_dist = far_dist,
        verbose = verbose
      )
    }
    else if (is.numeric(y)) {
      tsmessage(
        "Applying numeric set intersection, target weight = ",
        formatC(target_weight), " target neighbors = ", target_n_neighbors
      )
      target_nn <- find_nn(as.matrix(y), target_n_neighbors,
        method = "annoy",
        metric = "euclidean",
        n_trees = n_trees,
        n_threads = n_threads, grain_size = grain_size,
        search_k = search_k, verbose = FALSE
      )

      target_graph <- fuzzy_simplicial_set(
        nn = target_nn,
        set_op_mix_ratio = 1.0,
        local_connectivity = 1.0,
        bandwidth = 1.0,
        verbose = FALSE
      )

      V <- general_simplicial_set_intersection(
        V, target_graph, target_weight
      )

      V <- reset_local_connectivity(Matrix::drop0(V))
    }
    else {
      stop("y must be factors or numeric")
    }
  }


  if (methods::is(init, "matrix")) {
    if (nrow(init) != n_vertices || ncol(init) != n_components) {
      stop("init matrix does not match necessary configuration for X")
    }
    embedding <- init
  }
  else {
    init <- match.arg(tolower(init), c(
      "spectral", "random", "lvrandom", "normlaplacian",
      "laplacian", "spca", "pca"
    ))
    embedding <- switch(init,
      spectral = spectral_init(V, ndim = n_components, verbose = verbose),
      random = rand_init(n_vertices, n_components, verbose = verbose),
      lvrandom = rand_init_lv(n_vertices, n_components, verbose = verbose),
      normlaplacian = normalized_laplacian_init(V,
        ndim = n_components,
        verbose = verbose
      ),
      laplacian = laplacian_eigenmap(V, ndim = n_components, verbose = verbose),
      spca = scaled_pca(X, ndim = n_components, verbose = verbose),
      pca = pca_init(X, ndim = n_components, verbose = verbose),
      stop("Unknown initialization method: '", init, "'")
    )
  }


  if (is.null(n_epochs) || n_epochs <= 0) {
    if (method == "largevis") {
      n_epochs <- lvish_epochs(n_vertices, V)
    }
    else {
      if (n_vertices <= 10000) {
        n_epochs <- 500
      }
      else {
        n_epochs <- 200
      }
    }
  }

  V@x[V@x < max(V@x) / n_epochs] <- 0
  V <- Matrix::drop0(V)
  epochs_per_sample <- make_epochs_per_sample(V@x, n_epochs)

  positive_head <- V@i
  positive_tail <- Matrix::which(V != 0, arr.ind = TRUE)[, 2] - 1

  tsmessage(
    "Commencing optimization for ", n_epochs, " epochs, with ",
    length(positive_head), " positive edges",
    pluralize("thread", n_threads, " using")
  )

  parallelize <- n_threads > 0
  if (tolower(method) == "umap") {
    embedding <- optimize_layout_umap(
      head_embedding = embedding,
      tail_embedding = embedding,
      positive_head = positive_head,
      positive_tail = positive_tail,
      n_epochs = n_epochs,
      n_vertices = n_vertices,
      epochs_per_sample = epochs_per_sample,
      a = a, b = b, gamma = gamma,
      initial_alpha = alpha, negative_sample_rate,
      seed = get_seed(),
      approx_pow = approx_pow,
      parallelize = parallelize,
      grain_size = grain_size,
      move_other = TRUE,
      verbose = verbose
    )
  }
  else if (method == "tumap") {
    embedding <- optimize_layout_tumap(embedding,
      tail_embedding = embedding,
      positive_head = positive_head,
      positive_tail = positive_tail,
      n_epochs = n_epochs,
      n_vertices, epochs_per_sample,
      initial_alpha = alpha,
      negative_sample_rate = negative_sample_rate,
      seed = get_seed(),
      parallelize = parallelize,
      grain_size = grain_size,
      move_other = TRUE,
      verbose = verbose
    )
  }
  else {
    embedding <- optimize_layout_largevis(embedding,
      tail_embedding = embedding,
      positive_head = positive_head,
      positive_tail = positive_tail,
      n_epochs = n_epochs,
      n_vertices, epochs_per_sample,
      gamma = gamma,
      initial_alpha = alpha,
      negative_sample_rate = negative_sample_rate,
      seed = get_seed(),
      parallelize = parallelize,
      grain_size = grain_size,
      move_other = TRUE,
      verbose = verbose
    )
  }

  gc()
  # Center the points before returning
  embedding <- scale(embedding, center = TRUE, scale = FALSE)
  tsmessage("Optimization finished")

  if (ret_model) {
    list(
      scale_info = attr_to_scale_info(X),
      nn_index = nn$index,
      n_neighbors = n_neighbors,
      search_k = search_k,
      local_connectivity = local_connectivity,
      embedding = embedding,
      n_epochs = n_epochs,
      alpha = alpha,
      negative_sample_rate = negative_sample_rate,
      method = method,
      a = a,
      b = b,
      gamma = gamma,
      approx_pow = approx_pow
    )
  }
  else {
    embedding
  }
}



# Creates the number of epochs per sample for each weight
# weights are the non-zero input affinities (1-simplex)
# n_epoch the total number of epochs
# There is an inverse relationship between the weights and the return vector.
make_epochs_per_sample <- function(weights, n_epochs) {
  result <- rep(-1, length(weights))
  n_samples <- n_epochs * (weights / max(weights))
  result[n_samples > 0] <- n_epochs / n_samples[n_samples > 0]
  result
}

# Create the a/b parameters from spread and min_dist
find_ab_params <- function(spread = 1, min_dist = 0.001) {
  xv <- seq(from = 0, to = spread * 3, length.out = 300)
  yv <- rep(0, length(xv))
  yv[xv < min_dist] <- 1
  yv[xv >= min_dist] <- exp(-(xv[xv >= min_dist] - min_dist) / spread)
  result <- try({
    stats::nls(yv ~ 1 / (1 + a * xv^(2 * b)),
      start = list(a = 1, b = 1)
    )$m$getPars()
  }, silent = TRUE)
  if (class(result) == "try-error") {
    stop(
      "Can't find a, b for provided spread = ", spread,
      " min_dist = ", min_dist
    )
  }
  result
}

# The default number of edge samples used by LargeVis
lvish_samples <- function(n_vertices) {
  n_samples <- 0

  if (n_vertices < 10000) {
    n_samples <- 1000
  }
  else if (n_vertices < 1000000) {
    n_samples <- (n_vertices - 10000) * 9000 / (1000000 - 10000) + 1000
  }
  else {
    n_samples <- n_vertices / 100
  }

  round(n_samples * 1000000)
}

# Returns the number of epochs required to generate the default number of edge samples
# used in LargeVis
lvish_epochs <- function(n_vertices, V) {
  n_samples <- lvish_samples(n_vertices)
  round(n_samples * max(V) / sum(V))
}

# Scale X according to various strategies
scale_input <- function(X, scale_type, ret_model = FALSE, verbose = FALSE) {
  if (is.null(scale_type)) {
    scale_type <- "none"
  }
  else if (is.logical(scale_type)) {
    scale_type <- ifelse(scale_type, "scale", "none")
  }

  scale_type <- match.arg(
    tolower(scale_type),
    c("none", "scale", "range", "maxabs")
  )
  switch(scale_type,
    range = {
      tsmessage("Range scaling X")
      min_X <- min(X)
      X <- X - min_X

      max_X <- max(X)
      X <- X / max_X

      if (ret_model) {
        attr(X, "scaled:range:min") <- min_X
        attr(X, "scaled:range:max") <- max_X
      }
    },
    maxabs = {
      tsmessage("Normalizing by max-abs")
      X <- base::scale(X, scale = FALSE)
      max_abs <- max(abs(X))
      X <- X / max_abs

      if (ret_model) {
        attr(X, "scaled:maxabs") <- max_abs
      }
    },
    scale = {
      tsmessage("Scaling to zero mean and unit variance")

      varf <- function(x) {
        sum((x - sum(x) / length(x))^2)
      }
      non_zero_var_cols <- apply(X, 2, varf) >= .Machine$double.xmin

      if (length(non_zero_var_cols) == 0) {
        stop("Matrix has zero variance")
      }
      X <- X[, non_zero_var_cols]
      tsmessage("Kept ", ncol(X), " non-zero-variance columns")
      X <- base::scale(X, scale = TRUE)

      if (ret_model) {
        attr(X, "scaled:nzvcols") <- which(non_zero_var_cols)
      }
    }
  )
  X
}

attr_to_scale_info <- function(X) {
  Xattr <- attributes(X)
  Xattr <- Xattr[startsWith(names(Xattr), "scaled:")]
  if (length(Xattr) == 0) {
    Xattr <- NULL
  }
  Xattr
}

#' @useDynLib uwot, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
.onUnload <- function(libpath) {
  library.dynam.unload("uwot", libpath)
}
