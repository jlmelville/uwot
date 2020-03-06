#' Dimensionality Reduction with UMAP
#'
#' Carry out dimensionality reduction of a dataset using the Uniform Manifold
#' Approximation and Projection (UMAP) method (McInnes & Healy, 2018). Some of
#' the following help text is lifted verbatim from the Python reference
#' implementation at \url{https://github.com/lmcinnes/umap}.
#'
#' Note that the \code{grain_size} parameter no longer does anything and is
#' present to avoid break backwards compatibility only.
#'
#' @param X Input data. Can be a \code{\link{data.frame}}, \code{\link{matrix}},
#'   \code{\link[stats]{dist}} object or \code{\link[Matrix]{sparseMatrix}}. A
#'   sparse matrix is interpreted as a distance matrix and both implicit and
#'   explicit zero entries are ignored. Set zero distances you want to keep to
#'   an arbitrarily small non-zero value (e.g. \code{1e-10}). Matrix and data
#'   frames should contain one observation per row. Data frames will have any
#'   non-numeric columns removed, although factor columns will be used if
#'   explicitly included via \code{metric} (see the help for \code{metric} for
#'   details). Can be \code{NULL} if precomputed nearest neighbor data is passed
#'   to \code{nn_method}, and \code{init} is not \code{"spca"} or \code{"pca"}.
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
#'   \item \code{"hamming"}
#'   \item \code{"categorical"} (see below)
#' }
#' Only applies if \code{nn_method = "annoy"} (for \code{nn_method = "fnn"}, the
#' distance metric is always "euclidean").
#'
#' If \code{X} is a data frame or matrix, then multiple metrics can be
#' specified, by passing a list to this argument, where the name of each item in
#' the list is one of the metric names above. The value of each list item should
#' be a vector giving the names or integer ids of the columns to be included in
#' a calculation, e.g. \code{metric = list(euclidean = 1:4, manhattan = 5:10)}.
#'
#' Each metric calculation results in a separate fuzzy simplicial set, which are
#' intersected together to produce the final set. Metric names can be repeated.
#' Because non-numeric columns are removed from the data frame, it is safer to
#' use column names than integer ids.
#'
#' Factor columns can also be used by specifying the metric name
#' \code{"categorical"}. Factor columns are treated different from numeric
#' columns and although multiple factor columns can be specified in a vector,
#' each factor column specified is processed individually. If you specify
#' a non-factor column, it will be coerced to a factor.
#'
#' For a given data block, you may override the \code{pca} and \code{pca_center}
#' arguments for that block, by providing a list with one unnamed item
#' containing the column names or ids, and then any of the \code{pca} or
#' \code{pca_center} overrides as named items, e.g. \code{metric =
#' list(euclidean = 1:4, manhattan = list(5:10, pca_center = FALSE))}. This
#' exists to allow mixed binary and real-valued data to be included and to have
#' PCA applied to both, but with centering applied only to the real-valued data
#' (it is typical not to apply centering to binary data before PCA is applied).
#' @param n_epochs Number of epochs to use during the optimization of the
#'   embedded coordinates. By default, this value is set to \code{500} for
#'   datasets containing 10,000 vertices or less, and \code{200} otherwise.
#'   If \code{n_epochs = 0}, then coordinates determined by \code{"init"} will
#'   be returned.
#' @param scale Scaling to apply to \code{X} if it is a data frame or matrix:
#' \itemize{
#'   \item{\code{"none"} or \code{FALSE} or \code{NULL}} No scaling.
#'   \item{\code{"Z"} or \code{"scale"} or \code{TRUE}} Scale each column to
#'   zero mean and variance 1.
#'   \item{\code{"maxabs"}} Center each column to mean 0, then divide each
#'   element by the maximum absolute value over the entire matrix.
#'   \item{\code{"range"}} Range scale the entire matrix, so the smallest
#'   element is 0 and the largest is 1.
#'   \item{\code{"colrange"}} Scale each column in the range (0,1).
#' }
#' For UMAP, the default is \code{"none"}.
#' @param learning_rate Initial learning rate used in optimization of the
#'   coordinates.
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
#'     so the standard deviation is 1e-4, to give a distribution similar to that
#'     used in t-SNE. This is an alias for \code{init = "pca", init_sdev =
#'     1e-4}.
#'     \item \code{"agspectral"} An "approximate global" modification of
#'     \code{"spectral"} which all edges in the graph to a value of 1, and then
#'     sets a random number of edges (\code{negative_sample_rate} edges per
#'     vertex) to 0.1, to approximate the effect of non-local affinities.
#'     \item A matrix of initial coordinates.
#'   }
#'  For spectral initializations, (\code{"spectral"}, \code{"normlaplacian"},
#'  \code{"laplacian"}), if more than one connected component is identified,
#'  each connected component is initialized separately and the results are
#'  merged. If \code{verbose = TRUE} the number of connected components are
#'  logged to the console. The existence of multiple connected components
#'  implies that a global view of the data cannot be attained with this
#'  initialization. Either a PCA-based initialization or increasing the value of
#'  \code{n_neighbors} may be more appropriate.
#' @param init_sdev If non-\code{NULL}, scales each dimension of the initialized
#'   coordinates (including any user-supplied matrix) to this standard
#'   deviation. By default no scaling is carried out, except when \code{init =
#'   "spca"}, in which case the value is \code{0.0001}. Scaling the input may
#'   help if the unscaled versions result in initial coordinates with large
#'   inter-point distances or outliers. This usually results in small gradients
#'   during optimization and very little progress being made to the layout.
#'   Shrinking the initial embedding by rescaling can help under these
#'   circumstances. Scaling the result of \code{init = "pca"} is usually
#'   recommended and \code{init = "spca"} as an alias for \code{init = "pca",
#'   init_sdev = 1e-4} but for the spectral initializations the scaled versions
#'   usually aren't necessary unless you are using a large value of
#'   \code{n_neighbors} (e.g. \code{n_neighbors = 150} or higher).
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
#'   algorithm as similar to Laplacian Eigenmaps. Larger values induce more
#'   connectivity and a more global view of the data, smaller values concentrate
#'   more locally.
#' @param repulsion_strength Weighting applied to negative samples in low
#'   dimensional embedding optimization. Values higher than one will result in
#'   greater weight being given to negative samples.
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
#'   You may also pass precalculated nearest neighbor data to this argument. It
#'   must be a list consisting of two elements:
#'   \itemize{
#'     \item \code{"idx"}. A \code{n_vertices x n_neighbors} matrix
#'     containing the integer indexes of the nearest neighbors in \code{X}. Each
#'     vertex is considered to be its own nearest neighbor, i.e.
#'     \code{idx[, 1] == 1:n_vertices}.
#'     \item \code{"dist"}. A \code{n_vertices x n_neighbors} matrix
#'     containing the distances of the nearest neighbors.
#'   }
#'   Multiple nearest neighbor data (e.g. from two different precomputed
#'   metrics) can be passed by passing a list containing the nearest neighbor
#'   data lists as items.
#'   The \code{n_neighbors} parameter is ignored when using precomputed
#'   nearest neighbor data.
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
#' @param y Optional target data for supervised dimension reduction. Can be a
#' vector, matrix or data frame. Use the \code{target_metric} parameter to
#' specify the metrics to use, using the same syntax as \code{metric}. Usually
#' either a single numeric or factor column is used, but more complex formats
#' are possible. The following types are allowed:
#'   \itemize{
#'     \item Factor columns with the same length as \code{X}. \code{NA} is
#'     allowed for any observation with an unknown level, in which case
#'     UMAP operates as a form of semi-supervised learning. Each column is
#'     treated separately.
#'     \item Numeric data. \code{NA} is \emph{not} allowed in this case. Use the
#'     parameter \code{target_n_neighbors} to set the number of neighbors used
#'     with \code{y}. If unset, \code{n_neighbors} is used. Unlike factors,
#'     numeric columns are grouped into one block unless \code{target_metric}
#'     specifies otherwise. For example, if you wish columns \code{a} and
#'     \code{b} to be treated separately, specify
#'     \code{target_metric = list(euclidean = "a", euclidean = "b")}. Otherwise,
#'     the data will be effectively treated as a matrix with two columns.
#'     \item Nearest neighbor data, consisting of a list of two matrices,
#'     \code{idx} and \code{dist}. These represent the precalculated nearest
#'     neighbor indices and distances, respectively. This
#'     is the same format as that expected for precalculated data in
#'     \code{nn_method}. This format assumes that the underlying data was a
#'     numeric vector. Any user-supplied value of the \code{target_n_neighbors}
#'     parameter is ignored in this case, because the the number of columns in
#'     the matrices is used for the value. Multiple nearest neighbor data using
#'     different metrics can be supplied by passing a list of these lists.
#'   }
#' Unlike \code{X}, all factor columns included in \code{y} are automatically
#' used.
#' @param target_n_neighbors Number of nearest neighbors to use to construct the
#'   target simplicial set. Default value is \code{n_neighbors}. Applies only if
#'   \code{y} is non-\code{NULL} and \code{numeric}.
#' @param target_metric The metric used to measure distance for \code{y} if
#'   using supervised dimension reduction. Used only if \code{y} is numeric.
#' @param target_weight Weighting factor between data topology and target
#'   topology. A value of 0.0 weights entirely on data, a value of 1.0 weights
#'   entirely on target. The default of 0.5 balances the weighting equally
#'   between data and target. Only applies if \code{y} is non-\code{NULL}.
#' @param pca If set to a positive integer value, reduce data to this number of
#'   columns using PCA. Doesn't applied if the distance \code{metric} is
#'   \code{"hamming"}, or the dimensions of the data is larger than the
#'   number specified (i.e. number of rows and columns must be larger than the
#'   value of this parameter). If you have > 100 columns in a data frame or
#'   matrix, reducing the number of columns in this way may substantially
#'   increase the performance of the nearest neighbor search at the cost of a
#'   potential decrease in accuracy. In many t-SNE applications, a value of 50
#'   is recommended, although there's no guarantee that this is appropriate for
#'   all settings.
#' @param pca_center If \code{TRUE}, center the columns of \code{X} before
#'   carrying out PCA. For binary data, it's recommended to set this to
#'   \code{FALSE}.
#' @param pcg_rand If \code{TRUE}, use the PCG random number generator (O'Neill,
#'   2014) during optimization. Otherwise, use the faster (but probably less
#'   statistically good) Tausworthe "taus88" generator. The default is
#'   \code{TRUE}.
#' @param fast_sgd If \code{TRUE}, then the following combination of parameters
#'   is set: \code{pcg_rand = TRUE}, \code{n_sgd_threads = "auto"} and
#'   \code{approx_pow = TRUE}. The default is \code{FALSE}. Setting this to
#'   \code{TRUE} will speed up the stochastic optimization phase, but give a
#'   potentially less accurate embedding, and which will not be exactly
#'   reproducible even with a fixed seed. For visualization, \code{fast_sgd =
#'   TRUE} will give perfectly good results. For more generic dimensionality
#'   reduction, it's safer to leave \code{fast_sgd = FALSE}. If \code{fast_sgd =
#'   TRUE}, then user-supplied values of \code{pcg_rand}, \code{n_sgd_threads},
#'   and \code{approx_pow} are ignored.
#' @param ret_model If \code{TRUE}, then return extra data that can be used to
#'   add new data to an existing embedding via \code{\link{umap_transform}}. The
#'   embedded coordinates are returned as the list item \code{embedding}. If
#'   \code{FALSE}, just return the coordinates. This parameter can be used in
#'   conjunction with \code{ret_nn} and \code{ret_extra}. Note that some
#'   settings are incompatible with the production of a UMAP model: external
#'   neighbor data (passed via a list to \code{nn_method}), and factor columns
#'   that were included via the \code{metric} parameter. In the latter case, the
#'   model produced is based only on the numeric data. A transformation using
#'   new data is possible, but the factor columns in the new data are ignored.
#' @param ret_nn If \code{TRUE}, then in addition to the embedding, also return
#'   nearest neighbor data that can be used as input to \code{nn_method} to
#'   avoid the overhead of repeatedly calculating the nearest neighbors when
#'   manipulating unrelated parameters (e.g. \code{min_dist}, \code{n_epochs},
#'   \code{init}). See the "Value" section for the names of the list items. If
#'   \code{FALSE}, just return the coordinates. Note that the nearest neighbors
#'   could be sensitive to data scaling, so be wary of reusing nearest neighbor
#'   data if modifying the \code{scale} parameter. This parameter can be used in
#'   conjunction with \code{ret_model} and \code{ret_extra}.
#' @param ret_extra A vector indicating what extra data to return. May contain
#'   any combination of the following strings:
#'   \itemize{
#'     \item \code{"model"} Same as setting `ret_model = TRUE`.
#'     \item \code{"nn"} Same as setting `ret_nn = TRUE`.
#'     \item \code{"fgraph"} the high dimensional fuzzy graph (i.e. the fuzzy
#'       simplicial set of the merged local views of the input data). The graph
#'       is returned as a sparse symmetric N x N matrix of class
#'       \link[Matrix]{dgCMatrix-class}, where a non-zero entry (i, j) gives the
#'       membership strength of the edge connecting vertex i and vertex j. This
#'       can be considered analogous to the input probability (or similarity or
#'       affinity) used in t-SNE and LargeVis. Note that the graph is further
#'       sparsified by removing edges with sufficiently low membership strength
#'       that they would not be sampled by the probabilistic edge sampling
#'       employed for optimization and therefore the number of non-zero elements
#'       in the matrix is dependent on \code{n_epochs}. If you are only
#'       interested in the fuzzy input graph (e.g. for clustering), setting
#'       `n_epochs = 0` will avoid any further sparsifying.
#'   }
#' @param n_threads Number of threads to use (except during stochastic gradient
#'   descent). Default is half the number of concurrent threads supported by the
#'   system. For nearest neighbor search, only applies if 
#'   \code{nn_method = "annoy"}. If \code{n_threads > 1}, then the Annoy index 
#'   will be temporarily written to disk in the location determined by
#'   \code{\link[base]{tempfile}}.
#' @param n_sgd_threads Number of threads to use during stochastic gradient
#'   descent. If set to > 1, then results will not be reproducible, even if
#'   `set.seed` is called with a fixed seed before running. Set to
#'   \code{"auto"} go use the same value as \code{n_threads}.
#' @param grain_size The minimum amount of work to do on each thread. If this
#'   value is set high enough, then less than \code{n_threads} or
#'   \code{n_sgd_threads} will be used for processing, which might give a
#'   performance improvement if the overhead of thread management and context
#'   switching was outweighing the improvement due to concurrent processing.
#'   This should be left at default (\code{1}) and work will be spread evenly
#'   over all the threads specified.
#' @param tmpdir Temporary directory to store nearest neighbor indexes during
#'   nearest neighbor search. Default is \code{\link{tempdir}}. The index is
#'   only written to disk if \code{n_threads > 1} and
#'   \code{nn_method = "annoy"}; otherwise, this parameter is ignored.
#' @param verbose If \code{TRUE}, log details to the console.
#' @return A matrix of optimized coordinates, or:
#'   \itemize{
#'     \item if \code{ret_model = TRUE} (or \code{ret_extra} contains
#'     \code{"model"}), returns a list containing extra information that can be
#'     used to add new data to an existing embedding via
#'     \code{\link{umap_transform}}. In this case, the coordinates are available
#'     in the list item \code{embedding}.
#'     \item if \code{ret_nn = TRUE} (or \code{ret_extra} contains \code{"nn"}),
#'     returns the nearest neighbor data as a list called \code{nn}. This
#'     contains one list for each \code{metric} calculated, itself containing a
#'     matrix \code{idx} with the integer ids of the neighbors; and a matrix
#'     \code{dist} with the distances. The \code{nn} list (or a sub-list) can be
#'     used as input to the \code{nn_method} parameter.
#'     \item if \code{ret_extra} contains \code{"fgraph"} returns the high
#'     dimensional fuzzy graph as a sparse matrix called \code{fgraph}, of type
#'     \link[Matrix]{dgCMatrix-class}.
#'   }
#'   The returned list contains the combined data from any combination of
#'   specifying \code{ret_model}, \code{ret_nn} and \code{ret_extra}.
#' @examples
#' 
#' iris30 <- iris[c(1:10, 51:60, 101:110), ]
#' 
#' # Non-numeric columns are automatically removed so you can pass data frames
#' # directly in a lot of cases without pre-processing
#' iris_umap <- umap(iris30, n_neighbors = 5, learning_rate = 0.5, init = "random", n_epochs = 20)
#'
#' # Faster approximation to the gradient and return nearest neighbors
#' iris_umap <- umap(iris30, n_neighbors = 5, approx_pow = TRUE, ret_nn = TRUE, n_epochs = 20)
#'
#' # Can specify min_dist and spread parameters to control separation and size
#' # of clusters and reuse nearest neighbors for efficiency
#' nn <- iris_umap$nn
#' iris_umap <- umap(iris30, n_neighbors = 5, min_dist = 1, spread = 5, nn_method = nn, n_epochs = 20)
#'
#' # Supervised dimension reduction using the 'Species' factor column
#' iris_sumap <- umap(iris30, n_neighbors = 5, min_dist = 0.001, y = iris30$Species, 
#'                    target_weight = 0.5, n_epochs = 20)
#'
#' # Calculate Petal and Sepal neighbors separately (uses intersection of the resulting sets):
#' iris_umap <- umap(iris30, metric = list(
#'   "euclidean" = c("Sepal.Length", "Sepal.Width"),
#'   "euclidean" = c("Petal.Length", "Petal.Width")
#' ))
#'
#' @references
#' Belkin, M., & Niyogi, P. (2002).
#' Laplacian eigenmaps and spectral techniques for embedding and clustering.
#' In \emph{Advances in neural information processing systems}
#' (pp. 585-591).
#' \url{http://papers.nips.cc/paper/1961-laplacian-eigenmaps-and-spectral-techniques-for-embedding-and-clustering.pdf}
#'
#' McInnes, L., & Healy, J. (2018).
#' UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction
#' \emph{arXiv preprint} \emph{arXiv}:1802.03426.
#' \url{https://arxiv.org/abs/1802.03426}
#'
#' O’Neill, M. E. (2014).
#' \emph{PCG: A family of simple fast space-efficient statistically good
#' algorithms for random number generation}
#' (Report No. HMC-CS-2014-0905). Harvey Mudd College.
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
                 n_epochs = NULL, learning_rate = 1, scale = FALSE,
                 init = "spectral", init_sdev = NULL,
                 spread = 1, min_dist = 0.01,
                 set_op_mix_ratio = 1.0, local_connectivity = 1.0,
                 bandwidth = 1.0, repulsion_strength = 1.0,
                 negative_sample_rate = 5.0, a = NULL, b = NULL,
                 nn_method = NULL, n_trees = 50,
                 search_k = 2 * n_neighbors * n_trees,
                 approx_pow = FALSE,
                 y = NULL, target_n_neighbors = n_neighbors,
                 target_metric = "euclidean",
                 target_weight = 0.5,
                 pca = NULL, pca_center = TRUE,
                 pcg_rand = TRUE,
                 fast_sgd = FALSE,
                 ret_model = FALSE, ret_nn = FALSE, ret_extra = c(),
                 n_threads = NULL,
                 n_sgd_threads = 0,
                 grain_size = 1,
                 tmpdir = tempdir(),
                 verbose = getOption("verbose", TRUE)) {
  uwot(
    X = X, n_neighbors = n_neighbors, n_components = n_components,
    metric = metric, n_epochs = n_epochs, alpha = learning_rate, scale = scale,
    init = init, init_sdev = init_sdev,
    spread = spread, min_dist = min_dist,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity, bandwidth = bandwidth,
    gamma = repulsion_strength, negative_sample_rate = negative_sample_rate,
    a = a, b = b, nn_method = nn_method, n_trees = n_trees,
    search_k = search_k,
    method = "umap", approx_pow = approx_pow,
    n_threads = n_threads, n_sgd_threads = n_sgd_threads,
    grain_size = grain_size,
    y = y, target_n_neighbors = target_n_neighbors,
    target_weight = target_weight, target_metric = target_metric,
    pca = pca, pca_center = pca_center,
    pcg_rand = pcg_rand,
    fast_sgd = fast_sgd,
    ret_model = ret_model || "model" %in% ret_extra,
    ret_nn = ret_nn || "nn" %in% ret_extra,
    ret_fgraph = "fgraph" %in% ret_extra,
    tmpdir = tempdir(),
    verbose = verbose
  )
}

#' Dimensionality Reduction Using t-Distributed UMAP (t-UMAP)
#'
#' A faster (but less flexible) version of the UMAP gradient. For more detail on
#' UMAP, see the  \code{\link{umap}} function.
#'
#' By setting the UMAP curve parameters \code{a} and \code{b} to \code{1}, you
#' get back the Cauchy distribution as used in t-SNE and LargeVis. It also
#' results in a substantially simplified gradient expression. This can give
#' a speed improvement of around 50\%.
#' 
#' Note that the \code{grain_size} parameter no longer does anything and is
#' present to avoid break backwards compatibility only.
#'
#' @param X Input data. Can be a \code{\link{data.frame}}, \code{\link{matrix}},
#'   \code{\link[stats]{dist}} object or \code{\link[Matrix]{sparseMatrix}}. A
#'   sparse matrix is interpreted as a distance matrix and both implicit and
#'   explicit zero entries are ignored. Set zero distances you want to keep to
#'   an arbitrarily small non-zero value (e.g. \code{1e-10}). Matrix and data
#'   frames should contain one observation per row. Data frames will have any
#'   non-numeric columns removed, although factor columns will be used if
#'   explicitly included via \code{metric} (see the help for \code{metric} for
#'   details). Can be \code{NULL} if precomputed nearest neighbor data is passed
#'   to \code{nn_method}, and \code{init} is not \code{"spca"} or \code{"pca"}.
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
#'   \item \code{"hamming"}
#'   \item \code{"categorical"} (see below)
#' }
#' Only applies if \code{nn_method = "annoy"} (for \code{nn_method = "fnn"}, the
#' distance metric is always "euclidean").
#'
#' If \code{X} is a data frame or matrix, then multiple metrics can be
#' specified, by passing a list to this argument, where the name of each item in
#' the list is one of the metric names above. The value of each list item should
#' be a vector giving the names or integer ids of the columns to be included in
#' a calculation, e.g. \code{metric = list(euclidean = 1:4, manhattan = 5:10)}.
#'
#' Each metric calculation results in a separate fuzzy simplicial set, which are
#' intersected together to produce the final set. Metric names can be repeated.
#' Because non-numeric columns are removed from the data frame, it is safer to
#' use column names than integer ids.
#'
#' Factor columns can also be used by specifying the metric name
#' \code{"categorical"}. Factor columns are treated different from numeric
#' columns and although multiple factor columns can be specified in a vector,
#' each factor column specified is processed individually. If you specify
#' a non-factor column, it will be coerced to a factor.
#'
#' For a given data block, you may override the \code{pca} and \code{pca_center}
#' arguments for that block, by providing a list with one unnamed item
#' containing the column names or ids, and then any of the \code{pca} or
#' \code{pca_center} overrides as named items, e.g. \code{metric =
#' list(euclidean = 1:4, manhattan = list(5:10, pca_center = FALSE))}. This
#' exists to allow mixed binary and real-valued data to be included and to have
#' PCA applied to both, but with centering applied only to the real-valued data
#' (it is typical not to apply centering to binary data before PCA is applied).
#' @param n_epochs Number of epochs to use during the optimization of the
#'   embedded coordinates. By default, this value is set to \code{500} for
#'   datasets containing 10,000 vertices or less, and \code{200} otherwise.
#'   If \code{n_epochs = 0}, then coordinates determined by \code{"init"} will
#'   be returned.
#' @param learning_rate Initial learning rate used in optimization of the
#'   coordinates.
#' @param scale Scaling to apply to \code{X} if it is a data frame or matrix:
#' \itemize{
#'   \item{\code{"none"} or \code{FALSE} or \code{NULL}} No scaling.
#'   \item{\code{"Z"} or \code{"scale"} or \code{TRUE}} Scale each column to
#'   zero mean and variance 1.
#'   \item{\code{"maxabs"}} Center each column to mean 0, then divide each
#'   element by the maximum absolute value over the entire matrix.
#'   \item{\code{"range"}} Range scale the entire matrix, so the smallest
#'   element is 0 and the largest is 1.
#'   \item{\code{"colrange"}} Scale each column in the range (0,1).
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
#'     so the standard deviation is 1e-4, to give a distribution similar to that
#'     used in t-SNE. This is an alias for \code{init = "pca", init_sdev =
#'     1e-4}.
#'     \item \code{"agspectral"} An "approximate global" modification of
#'     \code{"spectral"} which all edges in the graph to a value of 1, and then
#'     sets a random number of edges (\code{negative_sample_rate} edges per
#'     vertex) to 0.1, to approximate the effect of non-local affinities.
#'     \item A matrix of initial coordinates.
#'   }
#'  For spectral initializations, (\code{"spectral"}, \code{"normlaplacian"},
#'  \code{"laplacian"}), if more than one connected component is identified,
#'  each connected component is initialized separately and the results are
#'  merged. If \code{verbose = TRUE} the number of connected components are
#'  logged to the console. The existence of multiple connected components
#'  implies that a global view of the data cannot be attained with this
#'  initialization. Either a PCA-based initialization or increasing the value of
#'  \code{n_neighbors} may be more appropriate.
#' @param init_sdev If non-\code{NULL}, scales each dimension of the initialized
#'   coordinates (including any user-supplied matrix) to this standard
#'   deviation. By default no scaling is carried out, except when \code{init =
#'   "spca"}, in which case the value is \code{0.0001}. Scaling the input may
#'   help if the unscaled versions result in initial coordinates with large
#'   inter-point distances or outliers. This usually results in small gradients
#'   during optimization and very little progress being made to the layout.
#'   Shrinking the initial embedding by rescaling can help under these
#'   circumstances. Scaling the result of \code{init = "pca"} is usually
#'   recommended and \code{init = "spca"} as an alias for \code{init = "pca",
#'   init_sdev = 1e-4} but for the spectral initializations the scaled versions
#'   usually aren't necessary unless you are using a large value of
#'   \code{n_neighbors} (e.g. \code{n_neighbors = 150} or higher).
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
#'   algorithm as similar to Laplacian Eigenmaps. Larger values induce more
#'   connectivity and a more global view of the data, smaller values concentrate
#'   more locally.
#' @param repulsion_strength Weighting applied to negative samples in low
#'   dimensional embedding optimization. Values higher than one will result in
#'   greater weight being given to negative samples.
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
#'   You may also pass precalculated nearest neighbor data to this argument. It
#'   must be a list consisting of two elements:
#'   \itemize{
#'     \item \code{"idx"}. A \code{n_vertices x n_neighbors} matrix
#'     containing the integer indexes of the nearest neighbors in \code{X}. Each
#'     vertex is considered to be its own nearest neighbor, i.e.
#'     \code{idx[, 1] == 1:n_vertices}.
#'     \item \code{"dist"}. A \code{n_vertices x n_neighbors} matrix
#'     containing the distances of the nearest neighbors.
#'   }
#'   Multiple nearest neighbor data (e.g. from two different precomputed
#'   metrics) can be passed by passing a list containing the nearest neighbor
#'   data lists as items.
#'   The \code{n_neighbors} parameter is ignored when using precalculated
#'   nearest neighbor data.
#' @param n_trees Number of trees to build when constructing the nearest
#'   neighbor index. The more trees specified, the larger the index, but the
#'   better the results. With \code{search_k}, determines the accuracy of the
#'   Annoy nearest neighbor search. Only used if the \code{nn_method} is
#'   \code{"annoy"}. Sensible values are between \code{10} to \code{100}.
#' @param search_k Number of nodes to search during the neighbor retrieval. The
#'   larger k, the more the accurate results, but the longer the search takes.
#'   With \code{n_trees}, determines the accuracy of the Annoy nearest neighbor
#'   search. Only used if the \code{nn_method} is \code{"annoy"}.
#' @param y Optional target data for supervised dimension reduction. Can be a
#' vector, matrix or data frame. Use the \code{target_metric} parameter to
#' specify the metrics to use, using the same syntax as \code{metric}. Usually
#' either a single numeric or factor column is used, but more complex formats
#' are possible. The following types are allowed:
#'   \itemize{
#'     \item Factor columns with the same length as \code{X}. \code{NA} is
#'     allowed for any observation with an unknown level, in which case
#'     UMAP operates as a form of semi-supervised learning. Each column is
#'     treated separately.
#'     \item Numeric data. \code{NA} is \emph{not} allowed in this case. Use the
#'     parameter \code{target_n_neighbors} to set the number of neighbors used
#'     with \code{y}. If unset, \code{n_neighbors} is used. Unlike factors,
#'     numeric columns are grouped into one block unless \code{target_metric}
#'     specifies otherwise. For example, if you wish columns \code{a} and
#'     \code{b} to be treated separately, specify
#'     \code{target_metric = list(euclidean = "a", euclidean = "b")}. Otherwise,
#'     the data will be effectively treated as a matrix with two columns.
#'     \item Nearest neighbor data, consisting of a list of two matrices,
#'     \code{idx} and \code{dist}. These represent the precalculated nearest
#'     neighbor indices and distances, respectively. This
#'     is the same format as that expected for precalculated data in
#'     \code{nn_method}. This format assumes that the underlying data was a
#'     numeric vector. Any user-supplied value of the \code{target_n_neighbors}
#'     parameter is ignored in this case, because the the number of columns in
#'     the matrices is used for the value. Multiple nearest neighbor data using
#'     different metrics can be supplied by passing a list of these lists.
#'   }
#' Unlike \code{X}, all factor columns included in \code{y} are automatically
#' used.
#' @param target_n_neighbors Number of nearest neighbors to use to construct the
#'   target simplicial set. Default value is \code{n_neighbors}. Applies only if
#'   \code{y} is non-\code{NULL} and \code{numeric}.
#' @param target_metric The metric used to measure distance for \code{y} if
#'   using supervised dimension reduction. Used only if \code{y} is numeric.
#' @param target_weight Weighting factor between data topology and target
#'   topology. A value of 0.0 weights entirely on data, a value of 1.0 weights
#'   entirely on target. The default of 0.5 balances the weighting equally
#'   between data and target. Only applies if \code{y} is non-\code{NULL}.
#' @param pca If set to a positive integer value, reduce data to this number of
#'   columns using PCA. Doesn't applied if the distance \code{metric} is
#'   \code{"hamming"}, or the dimensions of the data is larger than the
#'   number specified (i.e. number of rows and columns must be larger than the
#'   value of this parameter). If you have > 100 columns in a data frame or
#'   matrix, reducing the number of columns in this way may substantially
#'   increase the performance of the nearest neighbor search at the cost of a
#'   potential decrease in accuracy. In many t-SNE applications, a value of 50
#'   is recommended, although there's no guarantee that this is appropriate for
#'   all settings.
#' @param pca_center If \code{TRUE}, center the columns of \code{X} before
#'   carrying out PCA. For binary data, it's recommended to set this to
#'   \code{FALSE}.
#' @param pcg_rand If \code{TRUE}, use the PCG random number generator (O'Neill,
#'   2014) during optimization. Otherwise, use the faster (but probably less
#'   statistically good) Tausworthe "taus88" generator. The default is
#'   \code{TRUE}.
#' @param fast_sgd If \code{TRUE}, then the following combination of parameters
#'   is set: \code{pcg_rand = TRUE} and \code{n_sgd_threads = "auto"}. The
#'   default is \code{FALSE}. Setting this to \code{TRUE} will speed up the
#'   stochastic optimization phase, but give a potentially less accurate
#'   embedding, and which will not be exactly reproducible even with a fixed
#'   seed. For visualization, \code{fast_sgd = TRUE} will give perfectly good
#'   results. For more generic dimensionality reduction, it's safer to leave
#'   \code{fast_sgd = FALSE}. If \code{fast_sgd = TRUE}, then user-supplied
#'   values of \code{pcg_rand} and \code{n_sgd_threads}, are ignored.
#' @param ret_model If \code{TRUE}, then return extra data that can be used to
#'   add new data to an existing embedding via \code{\link{umap_transform}}. The
#'   embedded coordinates are returned as the list item \code{embedding}. If
#'   \code{FALSE}, just return the coordinates. This parameter can be used in
#'   conjunction with \code{ret_nn} and \code{ret_extra}. Note that some
#'   settings are incompatible with the production of a UMAP model: external
#'   neighbor data (passed via a list to \code{nn_method}), and factor columns
#'   that were included via the \code{metric} parameter. In the latter case, the
#'   model produced is based only on the numeric data. A transformation using
#'   new data is possible, but the factor columns in the new data are ignored.
#' @param ret_nn If \code{TRUE}, then in addition to the embedding, also return
#'   nearest neighbor data that can be used as input to \code{nn_method} to
#'   avoid the overhead of repeatedly calculating the nearest neighbors when
#'   manipulating unrelated parameters (e.g. \code{min_dist}, \code{n_epochs},
#'   \code{init}). See the "Value" section for the names of the list items. If
#'   \code{FALSE}, just return the coordinates. Note that the nearest neighbors
#'   could be sensitive to data scaling, so be wary of reusing nearest neighbor
#'   data if modifying the \code{scale} parameter. This parameter can be used in
#'   conjunction with \code{ret_model} and \code{ret_extra}.
#' @param ret_extra A vector indicating what extra data to return. May contain
#'   any combination of the following strings:
#'   \itemize{
#'     \item \code{"model"} Same as setting `ret_model = TRUE`.
#'     \item \code{"nn"} Same as setting `ret_nn = TRUE`.
#'     \item \code{"fgraph"} the high dimensional fuzzy graph (i.e. the fuzzy
#'       simplicial set of the merged local views of the input data). The graph
#'       is returned as a sparse symmetric N x N matrix of class
#'       \link[Matrix]{dgCMatrix-class}, where a non-zero entry (i, j) gives the
#'       membership strength of the edge connecting vertex i and vertex j. This
#'       can be considered analogous to the input probability (or similarity or
#'       affinity) used in t-SNE and LargeVis. Note that the graph is further
#'       sparsified by removing edges with sufficiently low membership strength
#'       that they would not be sampled by the probabilistic edge sampling
#'       employed for optimization and therefore the number of non-zero elements
#'       in the matrix is dependent on \code{n_epochs}. If you are only
#'       interested in the fuzzy input graph (e.g. for clustering), setting
#'       `n_epochs = 0` will avoid any further sparsifying.
#'   }
#' @param n_threads Number of threads to use (except during stochastic gradient
#'   descent). Default is half the number of concurrent threads supported by the
#'   system. For nearest neighbor search, only applies if 
#'   \code{nn_method = "annoy"}. If \code{n_threads > 1}, then the Annoy index 
#'   will be temporarily written to disk in the location determined by
#'   \code{\link[base]{tempfile}}.
#' @param n_sgd_threads Number of threads to use during stochastic gradient
#'   descent. If set to > 1, then results will not be reproducible, even if
#'   `set.seed` is called with a fixed seed before running. Set to
#'   \code{"auto"} go use the same value as \code{n_threads}.
#' @param grain_size The minimum amount of work to do on each thread. If this
#'   value is set high enough, then less than \code{n_threads} or
#'   \code{n_sgd_threads} will be used for processing, which might give a
#'   performance improvement if the overhead of thread management and context
#'   switching was outweighing the improvement due to concurrent processing.
#'   This should be left at default (\code{1}) and work will be spread evenly
#'   over all the threads specified.
#' @param tmpdir Temporary directory to store nearest neighbor indexes during
#'   nearest neighbor search. Default is \code{\link{tempdir}}. The index is
#'   only written to disk if \code{n_threads > 1} and
#'   \code{nn_method = "annoy"}; otherwise, this parameter is ignored.
#' @param verbose If \code{TRUE}, log details to the console.
#' @return A matrix of optimized coordinates, or:
#'   \itemize{
#'     \item if \code{ret_model = TRUE} (or \code{ret_extra} contains
#'     \code{"model"}), returns a list containing extra information that can be
#'     used to add new data to an existing embedding via
#'     \code{\link{umap_transform}}. In this case, the coordinates are available
#'     in the list item \code{embedding}.
#'     \item if \code{ret_nn = TRUE} (or \code{ret_extra} contains \code{"nn"}),
#'     returns the nearest neighbor data as a list called \code{nn}. This
#'     contains one list for each \code{metric} calculated, itself containing a
#'     matrix \code{idx} with the integer ids of the neighbors; and a matrix
#'     \code{dist} with the distances. The \code{nn} list (or a sub-list) can be
#'     used as input to the \code{nn_method} parameter.
#'     \item if \code{ret_extra} contains \code{"fgraph"} returns the high
#'     dimensional fuzzy graph as a sparse matrix called \code{fgraph}, of type
#'     \link[Matrix]{dgCMatrix-class}.
#'   }
#'   The returned list contains the combined data from any combination of
#'   specifying \code{ret_model}, \code{ret_nn} and \code{ret_extra}.
#' @examples
#' iris_tumap <- tumap(iris, n_neighbors = 50, learning_rate = 0.5)
#' @export
tumap <- function(X, n_neighbors = 15, n_components = 2, metric = "euclidean",
                  n_epochs = NULL,
                  learning_rate = 1, scale = FALSE,
                  init = "spectral", init_sdev = NULL,
                  set_op_mix_ratio = 1.0, local_connectivity = 1.0,
                  bandwidth = 1.0, repulsion_strength = 1.0,
                  negative_sample_rate = 5.0,
                  nn_method = NULL, n_trees = 50,
                  search_k = 2 * n_neighbors * n_trees,
                  n_threads = NULL,
                  n_sgd_threads = 0,
                  grain_size = 1,
                  y = NULL, target_n_neighbors = n_neighbors,
                  target_metric = "euclidean",
                  target_weight = 0.5,
                  pca = NULL, pca_center = TRUE,
                  pcg_rand = TRUE,
                  fast_sgd = FALSE,
                  ret_model = FALSE, ret_nn = FALSE, ret_extra = c(),
                  tmpdir = tempdir(),
                  verbose = getOption("verbose", TRUE)) {
  uwot(
    X = X, n_neighbors = n_neighbors, n_components = n_components,
    metric = metric,
    n_epochs = n_epochs, alpha = learning_rate, scale = scale,
    init = init, init_sdev = init_sdev,
    spread = NULL, min_dist = NULL, set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity, bandwidth = bandwidth,
    gamma = repulsion_strength, negative_sample_rate = negative_sample_rate,
    a = NULL, b = NULL, nn_method = nn_method, n_trees = n_trees,
    search_k = search_k,
    method = "tumap",
    n_threads = n_threads, n_sgd_threads = n_sgd_threads,
    grain_size = grain_size,
    y = y, target_n_neighbors = target_n_neighbors,
    target_weight = target_weight, target_metric = target_metric,
    pca = pca, pca_center = pca_center,
    pcg_rand = pcg_rand,
    fast_sgd = fast_sgd,
    ret_model = ret_model || "model" %in% ret_extra,
    ret_nn = ret_nn || "nn" %in% ret_extra,
    ret_fgraph = "fgraph" %in% ret_extra,
    tmpdir = tmpdir,
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
#' Note that the \code{grain_size} parameter no longer does anything and is
#' present to avoid break backwards compatibility only.
#'
#' @param X Input data. Can be a \code{\link{data.frame}}, \code{\link{matrix}},
#'   \code{\link[stats]{dist}} object or \code{\link[Matrix]{sparseMatrix}}. A
#'   sparse matrix is interpreted as a distance matrix and both implicit and
#'   explicit zero entries are ignored. Set zero distances you want to keep to
#'   an arbitrarily small non-zero value (e.g. \code{1e-10}). Matrix and data
#'   frames should contain one observation per row. Data frames will have any
#'   non-numeric columns removed, although factor columns will be used if
#'   explicitly included via \code{metric} (see the help for \code{metric} for
#'   details). Can be \code{NULL} if precomputed nearest neighbor data is passed
#'   to \code{nn_method}, and \code{init} is not \code{"spca"} or \code{"pca"}.
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
#'   \item \code{"hamming"}
#'   \item \code{"categorical"} (see below)
#' }
#' Only applies if \code{nn_method = "annoy"} (for \code{nn_method = "fnn"}, the
#' distance metric is always "euclidean").
#'
#' If \code{X} is a data frame or matrix, then multiple metrics can be
#' specified, by passing a list to this argument, where the name of each item in
#' the list is one of the metric names above. The value of each list item should
#' be a vector giving the names or integer ids of the columns to be included in
#' a calculation, e.g. \code{metric = list(euclidean = 1:4, manhattan = 5:10)}.
#'
#' Each metric calculation results in a separate fuzzy simplicial set, which are
#' intersected together to produce the final set. Metric names can be repeated.
#' Because non-numeric columns are removed from the data frame, it is safer to
#' use column names than integer ids.
#'
#' Factor columns can also be used by specifying the metric name
#' \code{"categorical"}. Factor columns are treated different from numeric
#' columns and although multiple factor columns can be specified in a vector,
#' each factor column specified is processed individually. If you specify
#' a non-factor column, it will be coerced to a factor.
#'
#' For a given data block, you may override the \code{pca} and \code{pca_center}
#' arguments for that block, by providing a list with one unnamed item
#' containing the column names or ids, and then any of the \code{pca} or
#' \code{pca_center} overrides as named items, e.g. \code{metric =
#' list(euclidean = 1:4, manhattan = list(5:10, pca_center = FALSE))}. This
#' exists to allow mixed binary and real-valued data to be included and to have
#' PCA applied to both, but with centering applied only to the real-valued data
#' (it is typical not to apply centering to binary data before PCA is applied).
#' @param n_epochs Number of epochs to use during the optimization of the
#'   embedded coordinates. The default is calculate the number of epochs
#'   dynamically based on dataset size, to give the same number of edge samples
#'   as the LargeVis defaults. This is usually substantially larger than the
#'   UMAP defaults. If \code{n_epochs = 0}, then coordinates determined by
#'   \code{"init"} will be returned.
#' @param learning_rate Initial learning rate used in optimization of the
#'   coordinates.
#' @param scale Scaling to apply to \code{X} if it is a data frame or matrix:
#' \itemize{
#'   \item{\code{"none"} or \code{FALSE} or \code{NULL}} No scaling.
#'   \item{\code{"Z"} or \code{"scale"} or \code{TRUE}} Scale each column to
#'   zero mean and variance 1.
#'   \item{\code{"maxabs"}} Center each column to mean 0, then divide each
#'   element by the maximum absolute value over the entire matrix.
#'   \item{\code{"range"}} Range scale the entire matrix, so the smallest
#'   element is 0 and the largest is 1.
#'   \item{\code{"colrange"}} Scale each column in the range (0,1).
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
#'     so the standard deviation is 1e-4, to give a distribution similar to that
#'     used in t-SNE and LargeVis. This is an alias for \code{init = "pca",
#'     init_sdev = 1e-4}.
#'     \item \code{"agspectral"} An "approximate global" modification of
#'     \code{"spectral"} which all edges in the graph to a value of 1, and then
#'     sets a random number of edges (\code{negative_sample_rate} edges per
#'     vertex) to 0.1, to approximate the effect of non-local affinities.
#'     \item A matrix of initial coordinates.
#'   }
#'  For spectral initializations, (\code{"spectral"}, \code{"normlaplacian"},
#'  \code{"laplacian"}), if more than one connected component is identified,
#'  each connected component is initialized separately and the results are
#'  merged. If \code{verbose = TRUE} the number of connected components are
#'  logged to the console. The existence of multiple connected components
#'  implies that a global view of the data cannot be attained with this
#'  initialization. Either a PCA-based initialization or increasing the value of
#'  \code{n_neighbors} may be more appropriate.
#' @param init_sdev If non-\code{NULL}, scales each dimension of the initialized
#'   coordinates (including any user-supplied matrix) to this standard
#'   deviation. By default no scaling is carried out, except when \code{init =
#'   "spca"}, in which case the value is \code{0.0001}. Scaling the input may
#'   help if the unscaled versions result in initial coordinates with large
#'   inter-point distances or outliers. This usually results in small gradients
#'   during optimization and very little progress being made to the layout.
#'   Shrinking the initial embedding by rescaling can help under these
#'   circumstances. Scaling the result of \code{init = "pca"} is usually
#'   recommended and \code{init = "spca"} as an alias for \code{init = "pca",
#'   init_sdev = 1e-4} but for the spectral initializations the scaled versions
#'   usually aren't necessary unless you are using a large value of
#'   \code{n_neighbors} (e.g. \code{n_neighbors = 150} or higher).
#' @param repulsion_strength Weighting applied to negative samples in low
#'   dimensional embedding optimization. Values higher than one will result in
#'   greater weight being given to negative samples.
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
#'   You may also pass precalculated nearest neighbor data to this argument. It
#'   must be a list consisting of two elements:
#'   \itemize{
#'     \item \code{"idx"}. A \code{n_vertices x n_neighbors} matrix
#'     containing the integer indexes of the nearest neighbors in \code{X}. Each
#'     vertex is considered to be its own nearest neighbor, i.e.
#'     \code{idx[, 1] == 1:n_vertices}.
#'     \item \code{"dist"}. A \code{n_vertices x n_neighbors} matrix
#'     containing the distances of the nearest neighbors.
#'   }
#'   Multiple nearest neighbor data (e.g. from two different precomputed
#'   metrics) can be passed by passing a list containing the nearest neighbor
#'   data lists as items.
#'   The \code{n_neighbors} parameter is ignored when using precomputed
#'   nearest neighbor data.
#' @param n_trees Number of trees to build when constructing the nearest
#'   neighbor index. The more trees specified, the larger the index, but the
#'   better the results. With \code{search_k}, determines the accuracy of the
#'   Annoy nearest neighbor search. Only used if the \code{nn_method} is
#'   \code{"annoy"}. Sensible values are between \code{10} to \code{100}.
#' @param search_k Number of nodes to search during the neighbor retrieval. The
#'   larger k, the more the accurate results, but the longer the search takes.
#'   With \code{n_trees}, determines the accuracy of the Annoy nearest neighbor
#'   search. Only used if the \code{nn_method} is \code{"annoy"}.
#' @param n_threads Number of threads to use (except during stochastic gradient
#'   descent). Default is half the number of concurrent threads supported by the
#'   system. For nearest neighbor search, only applies if 
#'   \code{nn_method = "annoy"}. If \code{n_threads > 1}, then the Annoy index 
#'   will be temporarily written to disk in the location determined by
#'   \code{\link[base]{tempfile}}.
#' @param n_sgd_threads Number of threads to use during stochastic gradient
#'   descent. If set to > 1, then results will not be reproducible, even if
#'   `set.seed` is called with a fixed seed before running. Set to
#'   \code{"auto"} go use the same value as \code{n_threads}.
#' @param grain_size The minimum amount of work to do on each thread. If this
#'   value is set high enough, then less than \code{n_threads} or
#'   \code{n_sgd_threads} will be used for processing, which might give a
#'   performance improvement if the overhead of thread management and context
#'   switching was outweighing the improvement due to concurrent processing.
#'   This should be left at default (\code{1}) and work will be spread evenly
#'   over all the threads specified.
#' @param kernel Type of kernel function to create input probabilities. Can be
#'   one of \code{"gauss"} (the default) or \code{"knn"}. \code{"gauss"} uses
#'   the usual Gaussian weighted similarities. \code{"knn"} assigns equal
#'   probabilities to every edge in the nearest neighbor graph, and zero
#'   otherwise, using \code{perplexity} nearest neighbors. The \code{n_neighbors}
#'   parameter is ignored in this case.
#' @param pca If set to a positive integer value, reduce data to this number of
#'   columns using PCA. Doesn't applied if the distance \code{metric} is
#'   \code{"hamming"}, or the dimensions of the data is larger than the
#'   number specified (i.e. number of rows and columns must be larger than the
#'   value of this parameter). If you have > 100 columns in a data frame or
#'   matrix, reducing the number of columns in this way may substantially
#'   increase the performance of the nearest neighbor search at the cost of a
#'   potential decrease in accuracy. In many t-SNE applications, a value of 50
#'   is recommended, although there's no guarantee that this is appropriate for
#'   all settings.
#' @param pca_center If \code{TRUE}, center the columns of \code{X} before
#'   carrying out PCA. For binary data, it's recommended to set this to
#'   \code{FALSE}.
#' @param pcg_rand If \code{TRUE}, use the PCG random number generator (O'Neill,
#'   2014) during optimization. Otherwise, use the faster (but probably less
#'   statistically good) Tausworthe "taus88" generator. The default is
#'   \code{TRUE}.
#' @param fast_sgd If \code{TRUE}, then the following combination of parameters
#'   is set: \code{pcg_rand = TRUE} and \code{n_sgd_threads = "auto"}. The
#'   default is \code{FALSE}. Setting this to \code{TRUE} will speed up the
#'   stochastic optimization phase, but give a potentially less accurate
#'   embedding, and which will not be exactly reproducible even with a fixed
#'   seed. For visualization, \code{fast_sgd = TRUE} will give perfectly good
#'   results. For more generic dimensionality reduction, it's safer to leave
#'   \code{fast_sgd = FALSE}. If \code{fast_sgd = TRUE}, then user-supplied
#'   values of \code{pcg_rand} and \code{n_sgd_threads}, are ignored.
#' @param ret_nn If \code{TRUE}, then in addition to the embedding, also return
#'   nearest neighbor data that can be used as input to \code{nn_method} to
#'   avoid the overhead of repeatedly calculating the nearest neighbors when
#'   manipulating unrelated parameters (e.g. \code{min_dist}, \code{n_epochs},
#'   \code{init}). See the "Value" section for the names of the list items. If
#'   \code{FALSE}, just return the coordinates. Note that the nearest neighbors
#'   could be sensitive to data scaling, so be wary of reusing nearest neighbor
#'   data if modifying the \code{scale} parameter.
#' @param ret_extra A vector indicating what extra data to return. May contain
#'   any combination of the following strings:
#'   \itemize{
#'     \item \code{"nn"} same as setting `ret_nn = TRUE`.
#'     \item \code{"P"} the high dimensional probability matrix. The graph
#'     is returned as a sparse symmetric N x N matrix of class
#'     \link[Matrix]{dgCMatrix-class}, where a non-zero entry (i, j) gives the
#'     input probability (or similarity or affinity) of the edge connecting
#'     vertex i and vertex j. Note that the graph is further sparsified by
#'     removing edges with sufficiently low membership strength that they would
#'     not be sampled by the probabilistic edge sampling employed for
#'     optimization and therefore the number of non-zero elements in the matrix
#'     is dependent on \code{n_epochs}. If you are only interested in the fuzzy
#'     input graph (e.g. for clustering), setting `n_epochs = 0` will avoid any
#'     further sparsifying.
#'   }
#' @param tmpdir Temporary directory to store nearest neighbor indexes during
#'   nearest neighbor search. Default is \code{\link{tempdir}}. The index is
#'   only written to disk if \code{n_threads > 1} and
#'   \code{nn_method = "annoy"}; otherwise, this parameter is ignored.
#' @param verbose If \code{TRUE}, log details to the console.
#' @return A matrix of optimized coordinates, or:
#'   \itemize{
#'     \item if \code{ret_nn = TRUE} (or \code{ret_extra} contains \code{"nn"}),
#'     returns the nearest neighbor data as a list called \code{nn}. This
#'     contains one list for each \code{metric} calculated, itself containing a
#'     matrix \code{idx} with the integer ids of the neighbors; and a matrix
#'     \code{dist} with the distances. The \code{nn} list (or a sub-list) can be
#'     used as input to the \code{nn_method} parameter.
#'     \item if \code{ret_extra} contains \code{"P"}, returns the high
#'     dimensional probability matrix as a sparse matrix called \code{P}, of
#'     type \link[Matrix]{dgCMatrix-class}.
#'   }
#'   The returned list contains the combined data from any combination of
#'   specifying \code{ret_nn} and \code{ret_extra}.
#' @references
#' Tang, J., Liu, J., Zhang, M., & Mei, Q. (2016, April).
#' Visualizing large-scale and high-dimensional data.
#' In \emph{Proceedings of the 25th International Conference on World Wide Web}
#' (pp. 287-297).
#' International World Wide Web Conferences Steering Committee.
#' \url{https://arxiv.org/abs/1602.00370}
#'
#' @examples
#' # Default number of epochs is much larger than for UMAP, assumes random
#' # initialization. Use perplexity rather than n_neighbors to control the size
#' # of the local neighborhood 20 epochs may be too small for a random 
#' # initialization
#' iris_lvish <- lvish(iris,
#'   perplexity = 50, learning_rate = 0.5,
#'   init = "random", n_epochs = 20
#' )
#' @export
lvish <- function(X, perplexity = 50, n_neighbors = perplexity * 3,
                  n_components = 2, metric = "euclidean", n_epochs = -1,
                  learning_rate = 1, scale = "maxabs",
                  init = "lvrandom", init_sdev = NULL,
                  repulsion_strength = 7,
                  negative_sample_rate = 5.0,
                  nn_method = NULL, n_trees = 50,
                  search_k = 2 * n_neighbors * n_trees,
                  n_threads = NULL,
                  n_sgd_threads = 0,
                  grain_size = 1,
                  kernel = "gauss",
                  pca = NULL, pca_center = TRUE,
                  pcg_rand = TRUE,
                  fast_sgd = FALSE,
                  ret_nn = FALSE, ret_extra = c(),
                  tmpdir = tempdir(),
                  verbose = getOption("verbose", TRUE)) {
  uwot(X,
    n_neighbors = n_neighbors, n_components = n_components,
    metric = metric,
    n_epochs = n_epochs, alpha = learning_rate, scale = scale,
    init = init, init_sdev = init_sdev,
    gamma = repulsion_strength, negative_sample_rate = negative_sample_rate,
    nn_method = nn_method,
    n_trees = n_trees, search_k = search_k,
    method = "largevis", perplexity = perplexity,
    pca = pca, pca_center = pca_center,
    n_threads = n_threads, n_sgd_threads = n_sgd_threads,
    grain_size = grain_size,
    kernel = kernel,
    ret_nn = ret_nn || "nn" %in% ret_extra,
    ret_fgraph = "P" %in% ret_extra,
    pcg_rand = pcg_rand,
    fast_sgd = fast_sgd,
    tmpdir = tmpdir,
    verbose = verbose
  )
}

# Function that does all the real work
uwot <- function(X, n_neighbors = 15, n_components = 2, metric = "euclidean",
                 n_epochs = NULL,
                 alpha = 1, scale = FALSE,
                 init = "spectral", init_sdev = NULL,
                 spread = 1, min_dist = 0.01,
                 set_op_mix_ratio = 1.0, local_connectivity = 1.0,
                 bandwidth = 1.0, gamma = 1.0,
                 negative_sample_rate = 5.0, a = NULL, b = NULL,
                 nn_method = NULL, n_trees = 50,
                 search_k = 2 * n_neighbors * n_trees,
                 method = "umap",
                 perplexity = 50, approx_pow = FALSE,
                 y = NULL, target_n_neighbors = n_neighbors,
                 target_metric = "euclidean",
                 target_weight = 0.5,
                 n_threads = NULL,
                 n_sgd_threads = 0,
                 grain_size = 1,
                 kernel = "gauss",
                 ret_model = FALSE, ret_nn = FALSE, ret_fgraph = FALSE,
                 pca = NULL, pca_center = TRUE,
                 pcg_rand = TRUE,
                 fast_sgd = FALSE,
                 tmpdir = tempdir(),
                 verbose = getOption("verbose", TRUE)) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  if (method == "umap" && (is.null(a) || is.null(b))) {
    ab_res <- find_ab_params(spread = spread, min_dist = min_dist)
    a <- ab_res[1]
    b <- ab_res[2]
    tsmessage("UMAP embedding parameters a = ", formatC(a), " b = ", formatC(b))
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

  if (!is.null(y) && is.numeric(y) && any(is.na(y))) {
    stop("numeric y cannot contain NA")
  }

  if (!is.numeric(n_components) || n_components < 1) {
    stop("'n_components' must be a positive integer")
  }

  if (!is.null(pca)) {
    if (!is.numeric(pca) || pca < 1) {
      stop("'pca' must be a positive integer")
    }
    if (pca < n_components) {
      stop("'pca' must be >= n_components")
    }
  }

  if (fast_sgd) {
    n_sgd_threads <- "auto"
    pcg_rand <- FALSE
    approx_pow <- TRUE
  }

  if (n_threads < 0) {
    stop("n_threads cannot be < 0")
  }
  if (n_threads %% 1 != 0) {
    n_threads <- round(n_threads)
    tsmessage("Non-integer 'n_threads' provided. Setting to ", n_threads)
  }
  if (n_sgd_threads == "auto") {
    n_sgd_threads <- n_threads
  }
  if (n_sgd_threads < 0) {
    stop("n_sgd_threads cannot be < 0")
  }
  if (n_sgd_threads %% 1 != 0) {
    n_sgd_threads <- round(n_sgd_threads)
    tsmessage("Non-integer 'n_sgd_threads' provided. Setting to ", n_sgd_threads)
  }

  # Store categorical columns to be used to generate the graph
  Xcat <- NULL
  # number of original columns in data frame (or matrix)
  # will be used only if using df or matrix and ret_model = TRUE
  norig_col <- NULL
  if (is.null(X)) {
    if (!is.list(nn_method)) {
      stop("If X is NULL, must provide NN data in nn_method")
    }
    if (is.character(init) && tolower(init) %in% c("spca", "pca")) {
      stop("init = 'pca' and 'spca' can't be used with X = NULL")
    }
    n_vertices <- x2nv(nn_method)
  }
  else if (methods::is(X, "dist")) {
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
    cat_ids <- NULL
    norig_col <- ncol(X)
    if (methods::is(X, "data.frame") || methods::is(X, "matrix")) {
      if (methods::is(X, "matrix")) {
        X <- data.frame(X)
      }
      cat_res <- find_categoricals(metric)
      metric <- cat_res$metrics
      cat_ids <- cat_res$categoricals
      # Convert categorical columns to factors if they aren't already
      if (!is.null(cat_ids)) {
        X[, cat_ids] <- lapply(X[, cat_ids, drop = FALSE], factor)
        Xcat <- X[, cat_ids, drop = FALSE]
      }

      indexes <- which(vapply(X, is.numeric, logical(1)))
      if (length(indexes) == 0) {
        stop("No numeric columns found")
      }
      X <- as.matrix(X[, indexes])
    }
    n_vertices <- nrow(X)
    tsmessage(
      "Read ", n_vertices, " rows and found ", ncol(X),
      " numeric columns",
      appendLF = is.null(cat_ids)
    )
    if (length(cat_ids) > 0) {
      tsmessage(" and ", pluralize("categorical column", length(cat_ids)),
        time_stamp = FALSE
      )
    }
    X <- scale_input(X,
      scale_type = scale, ret_model = ret_model,
      verbose = verbose
    )
  }

  if (method == "largevis" && kernel == "knn") {
    n_neighbors <- perplexity
  }

  if (n_neighbors > n_vertices) {
    # If nn_method is a list, we will determine n_neighbors later
    if (!is.list(nn_method)) {
      # Otherwise,for LargeVis, n_neighbors normally determined from perplexity
      # not an error to be too large
      if (method == "largevis") {
        tsmessage("Setting n_neighbors to ", n_vertices)
        n_neighbors <- n_vertices
      }
      else {
        stop("n_neighbors must be smaller than the dataset size")
      }
    }
  }

  if (!is.list(metric)) {
    metrics <- list(c())
    names(metrics) <- metric
  }
  else {
    metrics <- metric
  }

  # For typical case of numeric matrix X and not using hamming distance, save
  # PCA results here in case initialization uses PCA too
  pca_models <- NULL
  pca_shortcut <- FALSE
  if (!is.null(pca) && length(metric) == 1 && metric != "hamming" &&
    is.matrix(X) && ncol(X) > pca) {
    tsmessage("Reducing X column dimension to ", pca, " via PCA")
    pca_res <- pca_scores(X,
      ncol = pca, center = pca_center,
      ret_extra = ret_model, verbose = verbose
    )
    if (ret_model) {
      X <- pca_res$scores
      pca_models[["1"]] <- pca_res[c("center", "rotation")]
      pca_res <- NULL
    }
    else {
      X <- pca_res
    }
    pca_shortcut <- TRUE
  }

  d2sr <- data2set(X, Xcat, n_neighbors, metrics, nn_method,
    n_trees, search_k,
    method,
    set_op_mix_ratio, local_connectivity, bandwidth,
    perplexity, kernel,
    n_threads, grain_size,
    ret_model,
    pca = pca, pca_center = pca_center,
    n_vertices = n_vertices,
    tmpdir = tmpdir,
    verbose = verbose
  )

  V <- d2sr$V
  nns <- d2sr$nns
  if (is.null(pca_models)) {
    pca_models <- d2sr$pca_models
  }

  if (!is.null(y)) {
    tsmessage("Processing y data")

    if (!is.list(target_metric)) {
      target_metrics <- list(c())
      names(target_metrics) <- target_metric
    }
    else {
      target_metrics <- target_metric
    }

    ycat <- NULL
    ycat_ids <- NULL
    if (methods::is(y, "data.frame")) {
      ycat_res <- find_categoricals(target_metric)
      target_metric <- ycat_res$metrics
      ycat_ids <- ycat_res$categoricals
      if (!is.null(ycat_ids)) {
        ycat <- y[, ycat_ids, drop = FALSE]
      }
      else {
        ycindexes <- which(vapply(y, is.factor, logical(1)))
        if (length(ycindexes) > 0) {
          ycat <- (y[, ycindexes, drop = FALSE])
        }
      }

      yindexes <- which(vapply(y, is.numeric, logical(1)))

      if (length(yindexes) > 0) {
        y <- as.matrix(y[, yindexes])
      }
      else {
        y <- NULL
      }
    }
    else if (is.list(y)) {
      nn_method <- y
    }
    else if (is.numeric(y)) {
      y <- as.matrix(y)
    }
    else if (is.factor(y)) {
      ycat <- data.frame(y)
      y <- NULL
    }


    if (!is.null(y)) {
      yd2sr <- data2set(y, ycat, target_n_neighbors, target_metrics, nn_method,
        n_trees, search_k,
        method,
        set_op_mix_ratio = 1.0,
        local_connectivity = 1.0,
        bandwidth = 1.0,
        perplexity = perplexity, kernel = kernel,
        n_threads = n_threads, grain_size = grain_size,
        ret_model = FALSE,
        pca = pca,
        n_vertices = n_vertices,
        tmpdir = tmpdir,
        verbose = verbose
      )

      tsmessage(
        "Intersecting X and Y sets with target weight = ",
        formatC(target_weight)
      )
      V <- set_intersect(V, yd2sr$V, target_weight, reset = TRUE)
      yd2sr$V <- NULL
      yd2sr$nns <- NULL
    }
    else if (!is.null(ycat)) {
      V <- categorical_intersection_df(ycat, V,
        weight = target_weight,
        verbose = verbose
      )
    }
  }

  if (!(ret_model || ret_nn)) {
    nns <- NULL
    gc()
  }

  if (methods::is(init, "matrix")) {
    if (nrow(init) != n_vertices || ncol(init) != n_components) {
      stop("init matrix does not match necessary configuration for X: ", "should
           have dimensions (", n_vertices, ", ", n_components, ")")
    }
    tsmessage("Initializing from user-supplied matrix")
    embedding <- init
  }
  else if (!(methods::is(init, "character") && length(init) == 1)) {
    stop("init should be either a matrix or string describing the ",
         "initialization method")
  }
  else {
    init <- match.arg(tolower(init), c(
      "spectral", "random", "lvrandom", "normlaplacian",
      "laplacian", "spca", "pca", "inormlaplacian", "ispectral",
      "agspectral"
    ))

    if (init_is_spectral(init)) {
      connected <- connected_components(V)
      if (connected$n_components > 1) {
        tsmessage(
          "Found ", connected$n_components, " connected components, ",
          "falling back to 'spca' initialization with init_sdev = 1"
        )
        init <- "spca"
        init_sdev <- 1
      }
    }

    # Don't repeat PCA initialization if we've already done it once
    if (pca_shortcut && init %in% c("spca", "pca") && pca >= n_components) {
      embedding <- X[, 1:n_components]
      if (init == "spca") {
        tsmessage("Initializing from scaled PCA")
      }
      else {
        tsmessage("Initializing from PCA")
      }
    }
    else {
      embedding <- switch(init,
        spectral = spectral_init(V, ndim = n_components, verbose = verbose),
        random = rand_init(n_vertices, n_components, verbose = verbose),
        lvrandom = rand_init_lv(n_vertices, n_components, verbose = verbose),
        normlaplacian = normalized_laplacian_init(V,
          ndim = n_components,
          verbose = verbose
        ),
        laplacian = laplacian_eigenmap(V, ndim = n_components, verbose = verbose),
        # we handle scaling pca below
        spca = pca_init(X, ndim = n_components, verbose = verbose),
        pca = pca_init(X, ndim = n_components, verbose = verbose),
        ispectral = irlba_spectral_init(V, ndim = n_components, verbose = verbose),
        inormlaplacian = irlba_normalized_laplacian_init(V,
          ndim = n_components,
          verbose = verbose
        ),
        agspectral = agspectral_init(V,
          n_neg_nbrs = negative_sample_rate,
          ndim = n_components, verbose = verbose
        ),
        stop("Unknown initialization method: '", init, "'")
      )
    }
    if (!is.null(init_sdev) || init == "spca") {
      if (is.null(init_sdev)) {
        init_sdev <- 1e-4
      }
      embedding <- shrink_coords(embedding, init_sdev)
    }
  }

  if (is.null(n_epochs) || n_epochs < 0) {
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

  if (n_epochs > 0) {
    V@x[V@x < max(V@x) / n_epochs] <- 0
    V <- Matrix::drop0(V)

    epochs_per_sample <- make_epochs_per_sample(V@x, n_epochs)

    positive_head <- V@i
    positive_tail <- Matrix::which(V != 0, arr.ind = TRUE)[, 2] - 1

    tsmessage(
      "Commencing optimization for ", n_epochs, " epochs, with ",
      length(positive_head), " positive edges",
      pluralize("thread", n_sgd_threads, " using")
    )

    embedding <- t(embedding)
    if (tolower(method) == "umap") {
      embedding <- optimize_layout_umap(
        head_embedding = embedding,
        tail_embedding = NULL,
        positive_head = positive_head,
        positive_tail = positive_tail,
        n_epochs = n_epochs,
        n_vertices = n_vertices,
        epochs_per_sample = epochs_per_sample,
        a = a, b = b, gamma = gamma,
        initial_alpha = alpha, negative_sample_rate,
        approx_pow = approx_pow,
        pcg_rand = pcg_rand,
        n_threads = n_sgd_threads,
        grain_size = grain_size,
        move_other = TRUE,
        verbose = verbose
      )
    }
    else if (method == "tumap") {
      embedding <- optimize_layout_tumap(
        head_embedding = embedding,
        tail_embedding = NULL,
        positive_head = positive_head,
        positive_tail = positive_tail,
        n_epochs = n_epochs,
        n_vertices = n_vertices,
        epochs_per_sample = epochs_per_sample,
        initial_alpha = alpha,
        negative_sample_rate = negative_sample_rate,
        pcg_rand = pcg_rand,
        n_threads = n_sgd_threads,
        grain_size = grain_size,
        move_other = TRUE,
        verbose = verbose
      )
    }
    else {
      embedding <- optimize_layout_largevis(
        head_embedding = embedding,
        positive_head = positive_head,
        positive_tail = positive_tail,
        n_epochs = n_epochs,
        n_vertices = n_vertices,
        epochs_per_sample = epochs_per_sample,
        gamma = gamma,
        initial_alpha = alpha,
        negative_sample_rate = negative_sample_rate,
        pcg_rand = pcg_rand,
        n_threads = n_sgd_threads,
        grain_size = grain_size,
        verbose = verbose
      )
    }
    embedding <- t(embedding)
    gc()
    # Center the points before returning
    embedding <- scale(embedding, center = TRUE, scale = FALSE)
    tsmessage("Optimization finished")
  }

  if (ret_model || ret_nn || ret_fgraph) {
    nblocks <- length(nns)
    res <- list(embedding = embedding)
    if (ret_model) {
      res <- append(res, list(
        scale_info = attr_to_scale_info(X),
        n_neighbors = n_neighbors,
        search_k = search_k,
        local_connectivity = local_connectivity,
        n_epochs = n_epochs,
        alpha = alpha,
        negative_sample_rate = negative_sample_rate,
        method = method,
        a = a,
        b = b,
        gamma = gamma,
        approx_pow = approx_pow,
        metric = metrics,
        norig_col = norig_col,
        pcg_rand = pcg_rand
      ))
      if (nblocks > 1) {
        res$nn_index <- list()
        for (i in 1:nblocks) {
          res$nn_index[[i]] <- nns[[i]]$index
        }
      }
      else {
        res$nn_index <- nns[[1]]$index
        if (is.null(res$metric[[1]])) {
          # 31: Metric usually lists column indices or names, NULL means use all
          # of them, but for loading the NN index we need the number of 
          # columns explicitly (we don't have access to the column dimension of
          # the input data at load time)
          # To be sure of the dimensionality, fetch the first item from the 
          # index and see how many elements are in the returned vector.
          res$metric[[1]] <- list(ndim = length(res$nn_index$getItemsVector(0)))
        }
      }
      if (!is.null(pca_models)) {
        res$pca_models <- pca_models
      }
    }
    if (ret_nn) {
      res$nn <- list()
      for (i in 1:nblocks) {
        res$nn[[i]] <- list(idx = nns[[i]]$idx, dist = nns[[i]]$dist)
      }
      names(res$nn) <- names(nns)
    }
    if (ret_fgraph) {
      if (method == "largevis") {
        res$P <- V
      }
      else {
        res$fgraph <- V
      }
    }
  }
  else {
    res <- embedding
  }

  res
}

#' Save or Load a Model
#'
#' Functions to write a UMAP model to a file, and to restore.
#'
#' @param model a UMAP model create by \code{\link{umap}}.
#' @param file name of the file where the model is to be saved or read from.
#' @param unload if \code{TRUE}, unload all nearest neighbor indexes for the
#'   model. The \code{model} will no longer be valid for use in 
#'   \code{\link{umap_transform}} and the temporary working directory used 
#'   during model saving will be deleted. You will need to reload the model
#'   with `load_uwot` to use the model. If \code{FALSE}, then the model can
#'   be re-used without reloading, but you must manually unload the NN index 
#'   when you are finished using it if you want to delete the temporary working
#'   directory. To unload manually, use \code{\link{unload_uwot}}. The 
#'   absolute path of the working directory is found in the `mod_dir` item of 
#'   the return value. 
#' @param verbose if \code{TRUE}, log information to the console.
#' @return \code{model} with one extra item: `mod_dir`, which contains the path
#' to the working directory. If \code{unload = FALSE} then this directory still
#' exists after this function returns, and can be cleaned up with 
#' \code{\link{unload_uwot}}. If you don't care about cleaning up this 
#' directory, or \code{unload = TRUE}, then you can ignore the return value.
#' @examples
#' iris_train <- iris[c(1:10, 51:60), ]
#' iris_test <- iris[100:110, ]
#' 
#' # create model
#' model <- umap(iris_train, ret_model = TRUE, n_epochs = 20)
#' 
#' # save without unloading: this leaves behind a temporary working directory
#' model_file <- tempfile("iris_umap")
#' model <- save_uwot(model, file = model_file)
#' 
#' # The model can continue to be used
#' test_embedding <- umap_transform(iris_test, model)
#' 
#' # To manually unload the model from memory when finished and to clean up
#' # the working directory (this doesn't touch your model file)
#' unload_uwot(model)
#' 
#' # At this point, model cannot be used with umap_transform, this would fail:
#' # test_embedding2 <- umap_transform(iris_test, model)
#' 
#' # restore the model: this also creates a temporary working directory
#' model2 <- load_uwot(file = model_file)
#' test_embedding2 <- umap_transform(iris_test, model2)
#' 
#' # Unload and clean up the loaded model temp directory
#' unload_uwot(model2)
#' 
#' # clean up the model file
#' unlink(model_file)
#' 
#' # save with unloading: this deletes the temporary working directory but
#' # doesn't allow the model to be re-used
#' model3 <- umap(iris_train, ret_model = TRUE, n_epochs = 20)
#' model_file3 <- tempfile("iris_umap")
#' model3 <- save_uwot(model3, file = model_file3, unload = TRUE)
#' 
#' @seealso \code{\link{load_uwot}}, \code{\link{unload_uwot}}
#' @export
save_uwot <- function(model, file, unload = FALSE, verbose = FALSE) {
  if (!all_nn_indices_are_loaded(model)) {
    stop("cannot save: NN index is unloaded")
  }
  
  wd <- getwd()
  model_file <- abspath(file)
  if (file.exists(model_file)) {
    stop("model file ", model_file, " already exists")    
  }
  
  tryCatch(
    {
      # create directory to store files in
      mod_dir <- tempfile(pattern = "dir")
      tsmessage("Creating temp model dir ", mod_dir)
      dir.create(mod_dir)
      
      # create the tempdir/uwot subdirectory
      uwot_dir <- file.path(mod_dir, "uwot")
      tsmessage("Creating dir ", mod_dir)
      dir.create(uwot_dir)
      
      # save model file to tempdir/uwot/model
      model_tmpfname <- file.path(uwot_dir, "model")
      saveRDS(model, file = model_tmpfname)
      
      # save each nn index inside tempdir/uwot/model
      metrics <- names(model$metric)
      n_metrics <- length(metrics)
      for (i in 1:n_metrics) {
        nn_tmpfname <- file.path(uwot_dir, paste0("nn", i))
        if (n_metrics == 1) {
          model$nn_index$save(nn_tmpfname)
        }
        else {
          model$nn_index[[i]]$save(nn_tmpfname)
        }
      }
      
      # archive the files under the temp dir into the single target file
      # change directory so the archive only contains one directory
      tsmessage("Changing to ", mod_dir)
      setwd(mod_dir)
      tmp_model_file <- abspath(file)
      tsmessage("Creating ", tmp_model_file)
      utils::tar(tarfile = tmp_model_file, files = "uwot/")
    },
    finally = {
      setwd(wd)
      if (model_file != tmp_model_file) {
        tsmessage("Copying ", tmp_model_file, " to ", model_file)
        file.copy(from = tmp_model_file, to = model_file)
      }
      model$mod_dir <- mod_dir
      if (unload) {
        unload_uwot(model, cleanup = TRUE, verbose = verbose)
      }
    }
  )
  model
}

#' Save or Load a Model
#'
#' Functions to write a UMAP model to a file, and to restore.
#'
#' @param file name of the file where the model is to be saved or read from.
#' @param verbose if \code{TRUE}, log information to the console.
#' @return The model saved at \code{file}, for use with 
#' \code{\link{umap_transform}}. Additionally, it contains an extra item: 
#' `mod_dir`, which contains the path to the temporary working directory used 
#' during loading of the model. This directory cannot be removed until this
#' model has been unloaded by using \code{\link{unload_uwot}}.
#' @examples
#' iris_train <- iris[c(1:10, 51:60), ]
#' iris_test <- iris[100:110, ]
#' 
#' # create model
#' model <- umap(iris_train, ret_model = TRUE, n_epochs = 20)
#' 
#' # save without unloading: this leaves behind a temporary working directory
#' model_file <- tempfile("iris_umap")
#' model <- save_uwot(model, file = model_file)
#' 
#' # The model can continue to be used
#' test_embedding <- umap_transform(iris_test, model)
#' 
#' # To manually unload the model from memory when finished and to clean up
#' # the working directory (this doesn't touch your model file)
#' unload_uwot(model)
#' 
#' # At this point, model cannot be used with umap_transform, this would fail:
#' # test_embedding2 <- umap_transform(iris_test, model)
#' 
#' # restore the model: this also creates a temporary working directory
#' model2 <- load_uwot(file = model_file)
#' test_embedding2 <- umap_transform(iris_test, model2)
#' 
#' # Unload and clean up the loaded model temp directory
#' unload_uwot(model2)
#' 
#' # clean up the model file
#' unlink(model_file)
#' 
#' # save with unloading: this deletes the temporary working directory but
#' # doesn't allow the model to be re-used
#' model3 <- umap(iris_train, ret_model = TRUE, n_epochs = 20)
#' model_file3 <- tempfile("iris_umap")
#' model3 <- save_uwot(model3, file = model_file3, unload = TRUE)
#' 
#' @seealso \code{\link{save_uwot}}, \code{\link{unload_uwot}}
#' @export
load_uwot <- function(file, verbose = FALSE) {
  # create directory to store files in
  mod_dir <- tempfile(pattern = "dir")
  tsmessage("Creating temp directory ", mod_dir)
  dir.create(mod_dir)
  
  utils::untar(abspath(file), exdir = mod_dir, verbose = verbose)
  
  model_fname <- file.path(mod_dir, "uwot/model")
  if (!file.exists(model_fname)) {
    stop("Can't find model in ", file)
  }
  model <- readRDS(file = model_fname)
  
  metrics <- names(model$metric)
  n_metrics <- length(metrics)
  
  for (i in 1:n_metrics) {
    nn_fname <- file.path(mod_dir, paste0("uwot/nn", i))
    if (!file.exists(nn_fname)) {
      stop("Can't find nearest neighbor index ", nn_fname, " in ", file)
    }
    metric <- metrics[[i]]
    # 31: need to specify the index dimensionality when creating the index
    if (is.list(model$metric[[i]])) {
      # in case where there is only one metric, the value is a one-item list
      # named 'ndim' giving the number of dimensions directly: all columns
      # are used in this metric
      ndim <- model$metric[[i]]$ndim
    }
    else {
      # otherwise, metric specifies the name or index used for each metric,
      # so the dimension is the number of them
      ndim = length(model$metric[[i]])
    }
    ann <- create_ann(metric, ndim = ndim)
    ann$load(nn_fname)
    if (n_metrics == 1) {
      model$nn_index <- ann
    }
    else {
      model$nn_index[[i]] <- ann
    }
  }
  model$mod_dir <- mod_dir
  
  model
}


#' Unload a Model
#'
#' Unloads the UMAP model. This prevents the model being used with 
#' \code{\link{umap_transform}}, but allows the temporary working directory
#' associated with saving or loading the model to be removed.
#'
#' @param model a UMAP model create by \code{\link{umap}}.
#' @param cleanup if \code{TRUE}, attempt to delete the temporary working
#'   directory that was used in either the save or load of the model.
#' @param verbose if \code{TRUE}, log information to the console.
#'
#' @examples
#' iris_train <- iris[c(1:10, 51:60), ]
#' iris_test <- iris[100:110, ]
#' 
#' # create model
#' model <- umap(iris_train, ret_model = TRUE, n_epochs = 20)
#' 
#' # save without unloading: this leaves behind a temporary working directory
#' model_file <- tempfile("iris_umap")
#' model <- save_uwot(model, file = model_file)
#' 
#' # The model can continue to be used
#' test_embedding <- umap_transform(iris_test, model)
#' 
#' # To manually unload the model from memory when finished and to clean up
#' # the working directory (this doesn't touch your model file)
#' unload_uwot(model)
#' 
#' # At this point, model cannot be used with umap_transform, this would fail:
#' # test_embedding2 <- umap_transform(iris_test, model)
#' 
#' # restore the model: this also creates a temporary working directory
#' model2 <- load_uwot(file = model_file)
#' test_embedding2 <- umap_transform(iris_test, model2)
#' 
#' # Unload and clean up the loaded model temp directory
#' unload_uwot(model2)
#' 
#' # clean up the model file
#' unlink(model_file)
#' 
#' # save with unloading: this deletes the temporary working directory but
#' # doesn't allow the model to be re-used
#' model3 <- umap(iris_train, ret_model = TRUE, n_epochs = 20)
#' model_file3 <- tempfile("iris_umap")
#' model3 <- save_uwot(model3, file = model_file3, unload = TRUE)
#'
#' @seealso \code{\link{save_uwot}}, \code{\link{load_uwot}}
#' @export
unload_uwot <- function(model, cleanup = TRUE, verbose = FALSE) {
  tsmessage("Unloading NN index: model will be invalid")
  metrics <- names(model$metric)
  n_metrics <- length(metrics)
  for (i in 1:n_metrics) {
    if (n_metrics == 1) {
      model$nn_index$unload()
    }
    else {
      model$nn_index[[i]]$unload()
    }
  }
  
  if (cleanup) {
    if (is.null(model$mod_dir)) {
      tsmessage("Model is missing temp dir location, can't clean up")
      return();
    }
    else {
      mod_dir <- model$mod_dir
      if (!file.exists(mod_dir)) {
        tsmessage("model temp dir location '", mod_dir, "' no longer exists")
        return()
      }
      tsmessage("Deleting temp model dir ", mod_dir)
      res <- unlink(mod_dir, recursive = TRUE)
      if (res != 0) {
        tsmessage("Unable to delete tempdir ", mod_dir)
      }
    }
  }
}

all_nn_indices_are_loaded <- function(model) {
  if (is.null(model$nn_index)) {
    stop("Invalid model: has no 'nn_index'")
  }
  if (is.list(model$nn_index)) {
    for (i in 1:length(model$nn_index)) {
      if (model$nn_index[[i]]$getNTrees() == 0) {
        return(FALSE)
      }
    }
  }
  else if (model$nn_index$getNTrees() == 0) {
    return(FALSE)
  }
  TRUE
}

abspath <- function(filename) {
  file.path(normalizePath(dirname(filename)), basename(filename))
}

# Half of whatever the C++ implementation thinks are the number of concurrent 
# threads supported, but at least 1
default_num_threads <- function() {
  max(1, hardware_concurrency() / 2)
}

# Get the number of vertices in X
x2nv <- function(X) {
  if (is.list(X)) {
    if (!is.null(X$idx)) {
      n_vertices <- x2nv(X$idx)
    }
    else {
      if (length(X) > 0) {
        n_vertices <- x2nv(X[[1]])
      }
      else {
        stop("Can't find n_vertices for list X")
      }
    }
  }
  else if (methods::is(X, "dist")) {
    n_vertices <- attr(X, "Size")
  }
  else if (methods::is(X, "sparseMatrix")) {
    n_vertices <- nrow(X)
  }
  else if (methods::is(X, "data.frame") || methods::is(X, "matrix")) {
    n_vertices <- nrow(X)
  }
  else if (is.numeric(X)) {
    n_vertices <- length(X)
  }
  else {
    stop("Can't find number of vertices for X of type '", class(X)[1], "'")
  }
  n_vertices
}

data2set <- function(X, Xcat, n_neighbors, metrics, nn_method,
                     n_trees, search_k,
                     method,
                     set_op_mix_ratio, local_connectivity, bandwidth,
                     perplexity, kernel,
                     n_threads, grain_size,
                     ret_model,
                     n_vertices = x2nv(X),
                     tmpdir = tempdir(),
                     pca = NULL, pca_center = TRUE,
                     verbose = FALSE) {
  V <- NULL
  nns <- list()
  nblocks <- length(metrics)

  # Check for precalculated NN data in nn_method
  if (is.list(nn_method)) {
    if (is.null(nn_method$idx)) {
      nblocks <- length(nn_method)
      if (nblocks == 0) {
        stop("Incorrect format for precalculated neighbor data")
      }
    }
    else {
      nblocks <- 1
      # wrap nn data in a list so data is always a list of lists
      nn_method <- list(nn_method)
    }
    metrics <- replicate(nblocks, NULL, simplify = FALSE)
    names(metrics) <- rep("precomputed", nblocks)
  }


  if (nblocks > 1) {
    tsmessage("Found ", nblocks, " blocks of data")
  }
  mnames <- tolower(names(metrics))
  if (is.null(nn_method)) {
    if (n_vertices < 4096 && !ret_model && all(mnames == "euclidean")) {
      tsmessage("Using FNN for neighbor search, n_neighbors = ", n_neighbors)
      nn_method <- "fnn"
    }
    else {
      tsmessage("Using Annoy for neighbor search, n_neighbors = ", n_neighbors)
      nn_method <- "annoy"
    }
  }

  pca_models <- list()
  for (i in 1:nblocks) {
    metric <- mnames[[i]]
    metric <- match.arg(metric, c(
      "euclidean", "cosine", "manhattan",
      "hamming", "precomputed"
    ))
    # Defaults for this block which can be overridden
    pca_i <- pca
    pca_center_i <- pca_center

    subset <- metrics[[i]]
    if (is.null(subset)) {
      Xsub <- X
    }
    else if (is.list(subset)) {
      # e.g. "euclidean" = list(1:10, pca_center = FALSE),
      lsres <- lsplit_unnamed(subset)
      if (is.null(lsres$unnamed)) {
        stop("Error: no subset provided for block ", i)
      }
      if (length(lsres$unnamed) != 1) {
        stop("Error: only one unnamed item should be provided for block ", i)
      }
      subset <- lsres$unnamed[[1]]

      # possible overrides
      if (!is.null(lsres$named)) {
        lsnamed <- lsres$named
        lsnames <- names(lsnamed)
        if (!is.null(lsnamed$pca_center)) {
          pca_center_i <- lsnamed$pca_center
        }
        # PCA argument can be NULL, so need to check if it was explicitly provided
        if ("pca" %in% lsnames) {
          pca_i <- lsnamed$pca
        }
      }
      Xsub <- X[, subset, drop = FALSE]
    }
    else {
      Xsub <- X[, subset, drop = FALSE]
    }

    if (!is.null(X) && is.matrix(X)) {
      block_size <- ncol(Xsub)
      if (block_size == 0) {
        stop("Block ", i, " has zero size")
      }
      if (nblocks > 1) {
        tsmessage(
          "Processing block ", i, " of ", nblocks,
          " with size ", block_size,
          " using metric '", metric, "'"
        )
      }
    }
    else {
      # X is NULL or dist or something like that
      if (nblocks > 1) {
        tsmessage(
          "Processing block ", i, " of ", nblocks,
          " using metric '", metric, "'"
        )
      }
    }

    if (!is.null(pca_i) && is.matrix(X) && metric != "hamming" &&
      ncol(X) > pca_i && nrow(X) > pca_i) {
      tsmessage("Reducing column dimension to ", pca_i, " via PCA")
      pca_res <- pca_scores(Xsub, pca_i,
        ret_extra = ret_model,
        center = pca_center_i,
        verbose = verbose
      )
      if (ret_model) {
        Xsub <- pca_res$scores
        pca_models[[as.character(i)]] <- pca_res[c("center", "rotation")]
        pca_res <- NULL
      }
      else {
        Xsub <- pca_res
      }
    }

    nn_sub <- nn_method
    # Extract this block of nn data from list of lists
    if (metric == "precomputed") {
      nn_sub <- nn_method[[i]]
      if (i == 1) {
        n_neighbors <- NULL
      }
      else {
        n_neighbors <- ncol(nn_method[[1]]$idx)
      }
    }

    x2set_res <- x2set(Xsub, n_neighbors, metric,
      nn_method = nn_sub,
      n_trees, search_k,
      method,
      set_op_mix_ratio, local_connectivity, bandwidth,
      perplexity, kernel,
      n_threads, grain_size,
      ret_model,
      n_vertices = n_vertices,
      tmpdir = tmpdir,
      verbose = verbose
    )
    Vblock <- x2set_res$V
    nn <- x2set_res$nn
    nns[[i]] <- nn
    names(nns)[[i]] <- metric
    n_neighbors <- ncol(nn$idx)
    if (is.null(V)) {
      V <- Vblock
    }
    else {
      V <- set_intersect(V, Vblock, weight = 0.5, reset = TRUE)
    }
  }

  if (!is.null(Xcat)) {
    V <- categorical_intersection_df(Xcat, V, weight = 0.5, verbose = verbose)
  }

  list(V = V, nns = nns, pca_models = pca_models)
}

x2nn <- function(X, n_neighbors, metric, nn_method,
                 n_trees, search_k,
                 tmpdir = tempdir(),
                 n_threads, grain_size,
                 ret_model,
                 n_vertices = x2nv(X),
                 verbose = FALSE) {
  if (is.list(nn_method)) {
    # on first iteration n_neighbors is NULL
    # on subsequent iterations ensure n_neighbors is consistent for all data
    validate_nn(nn_method, n_vertices, n_neighbors = n_neighbors)
    nn <- nn_method
  }
  else {
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
      n_trees = n_trees, search_k = search_k,
      tmpdir = tmpdir,
      n_threads = n_threads, grain_size = grain_size,
      ret_index = ret_model, verbose = verbose
    )
  }
  nn
}

validate_nn <- function(nn_method, n_vertices, n_neighbors = NULL) {
  if (!is.matrix(nn_method$idx)) {
    stop("Couldn't find precalculated 'idx' matrix")
  }
  if (nrow(nn_method$idx) != n_vertices) {
    stop(
      "Precalculated 'idx' matrix must have ", n_vertices,
      " rows, but found ", nrow(nn_method$idx)
    )
  }

  # set n_neighbors from these matrices if it hasn't been already set
  if (is.null(n_neighbors)) {
    n_neighbors <- ncol(nn_method$idx)
  }
  if (!is.matrix(nn_method$dist)) {
    stop("Couldn't find precalculated 'dist' matrix")
  }
  if (nrow(nn_method$idx) != n_vertices) {
    stop("Precalculated 'dist' matrix must have ", n_vertices, " rows, but
         found ", nrow(nn_method$dist))
  }
  if (ncol(nn_method$dist) != n_neighbors) {
    stop("Precalculated 'dist' matrix must have ", n_neighbors, " cols, but
         found ", ncol(nn_method$dist))
  }
}

nn2set <- function(method, nn,
                   set_op_mix_ratio, local_connectivity, bandwidth,
                   perplexity, kernel,
                   n_threads, grain_size,
                   verbose = FALSE) {
  if (method == "largevis") {
    n_vertices <- nrow(nn$dist)
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
}


x2set <- function(X, n_neighbors, metric, nn_method,
                  n_trees, search_k,
                  method,
                  set_op_mix_ratio, local_connectivity, bandwidth,
                  perplexity, kernel,
                  n_threads, grain_size,
                  ret_model,
                  n_vertices = x2nv(X),
                  tmpdir = tempdir(),
                  verbose = FALSE) {
  nn <- x2nn(X,
    n_neighbors = n_neighbors,
    metric = metric,
    nn_method = nn_method,
    n_trees = n_trees, search_k = search_k,
    tmpdir = tmpdir,
    n_threads = n_threads, grain_size = grain_size,
    ret_model = ret_model,
    n_vertices = n_vertices,
    verbose = verbose
  )

  if (any(is.infinite(nn$dist))) {
    stop("Infinite distances found in nearest neighbors")
  }
  gc()


  V <- nn2set(method, nn,
    set_op_mix_ratio, local_connectivity, bandwidth,
    perplexity, kernel,
    n_threads, grain_size,
    verbose = verbose
  )
  if (any(is.na(V))) {
    stop("Non-finite entries in the input matrix")
  }
  gc()

  list(
    nn = nn,
    V = V
  )
}

set_intersect <- function(A, B, weight = 0.5, reset = TRUE) {
  A <- general_simplicial_set_intersection(
    A, B, weight
  )

  # https://github.com/lmcinnes/umap/issues/58#issuecomment-437633658
  # For now always reset
  if (reset) {
    A <- reset_local_connectivity(Matrix::drop0(A))
  }
  A
}

categorical_intersection_df <- function(X, V, weight = 0.5, verbose = FALSE) {
  tsmessage(
    "Carrying out categorical intersection for ",
    pluralize("column", ncol(X))
  )
  for (i in 1:ncol(X)) {
    V <- categorical_intersection(X[, i], V,
      weight = weight,
      verbose = (verbose && i == 1)
    )
  }
  V
}

categorical_intersection <- function(x, V, weight, verbose = FALSE) {
  if (is.null(V)) {
    stop("V cannot be null for categorical intersection")
  }
  if (weight < 1.0) {
    far_dist <- 2.5 * (1.0 / (1.0 - weight))
  }
  else {
    far_dist <- 1.0e12
  }
  tsmessage(
    "Applying categorical set intersection, weight = ", formatC(weight),
    " far distance = ", formatC(far_dist)
  )

  V <- categorical_simplicial_set_intersection(V, x,
    far_dist = far_dist,
    verbose = verbose
  )
  V
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
  result <- try(
    {
      stats::nls(yv ~ 1 / (1 + a * xv^(2 * b)),
        start = list(a = 1, b = 1)
      )$m$getPars()
    },
    silent = TRUE
  )
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
  else if (tolower(scale_type) == "z") {
    scale_type <- "scale"
  }

  scale_type <- match.arg(
    tolower(scale_type),
    c("none", "scale", "range", "colrange", "maxabs")
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
    colrange = {
      tsmessage("Column range scaling X")
      min_X <- apply(X, 2, min)
      X <- sweep(X, 2, min_X)

      max_X <- apply(X, 2, max)
      X <- sweep(X, 2, max_X, `/`)

      if (ret_model) {
        attr(X, "scaled:colrange:min") <- min_X
        attr(X, "scaled:colrange:max") <- max_X
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
.onUnload <- function(libpath) {
  library.dynam.unload("uwot", libpath)
}
