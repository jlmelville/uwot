# Dimensionality Reduction with a LargeVis-like method

Carry out dimensionality reduction of a dataset using a method similar
to LargeVis (Tang et al., 2016).

## Usage

``` r
lvish(
  X,
  perplexity = 50,
  n_neighbors = perplexity * 3,
  n_components = 2,
  metric = "euclidean",
  n_epochs = -1,
  learning_rate = 1,
  scale = "maxabs",
  init = "lvrandom",
  init_sdev = NULL,
  repulsion_strength = 7,
  negative_sample_rate = 5,
  nn_method = NULL,
  n_trees = 50,
  search_k = 2 * n_neighbors * n_trees,
  n_threads = NULL,
  n_sgd_threads = 0,
  grain_size = 1,
  kernel = "gauss",
  pca = NULL,
  pca_center = TRUE,
  pcg_rand = TRUE,
  fast_sgd = FALSE,
  ret_nn = FALSE,
  ret_extra = c(),
  tmpdir = tempdir(),
  verbose = getOption("verbose", TRUE),
  batch = FALSE,
  opt_args = NULL,
  epoch_callback = NULL,
  pca_method = NULL,
  binary_edge_weights = FALSE,
  nn_args = list(),
  rng_type = NULL
)
```

## Arguments

- X:

  Input data. Can be a
  [`data.frame`](https://rdrr.io/r/base/data.frame.html),
  [`matrix`](https://rdrr.io/r/base/matrix.html),
  [`dist`](https://rdrr.io/r/stats/dist.html) object or
  [`sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html).
  Matrix and data frames should contain one observation per row. Data
  frames will have any non-numeric columns removed, although factor
  columns will be used if explicitly included via `metric` (see the help
  for `metric` for details). A sparse matrix is interpreted as a
  distance matrix, and is assumed to be symmetric, so you can also pass
  in an explicitly upper or lower triangular sparse matrix to save
  storage. There must be at least `n_neighbors` non-zero distances for
  each row. Both implicit and explicit zero entries are ignored. Set
  zero distances you want to keep to an arbitrarily small non-zero value
  (e.g. `1e-10`). `X` can also be `NULL` if pre-computed nearest
  neighbor data is passed to `nn_method`, and `init` is not `"spca"` or
  `"pca"`.

- perplexity:

  Controls the size of the local neighborhood used for manifold
  approximation. This is the analogous to `n_neighbors` in
  [`umap`](https://jlmelville.github.io/uwot/reference/umap.md). Change
  this, rather than `n_neighbors`.

- n_neighbors:

  The number of neighbors to use when calculating the `perplexity`.
  Usually set to three times the value of the `perplexity`. Must be at
  least as large as `perplexity`.

- n_components:

  The dimension of the space to embed into. This defaults to `2` to
  provide easy visualization, but can reasonably be set to any integer
  value in the range `2` to `100`.

- metric:

  Type of distance metric to use to find nearest neighbors. For
  `nn_method = "annoy"` this can be one of:

  - `"euclidean"` (the default)

  - `"cosine"`

  - `"manhattan"`

  - `"hamming"`

  - `"correlation"` (a distance based on the Pearson correlation)

  - `"categorical"` (see below)

  For `nn_method = "hnsw"` this can be one of:

  - `"euclidean"`

  - `"cosine"`

  - `"correlation"`

  If [rnndescent](https://cran.r-project.org/package=rnndescent) is
  installed and `nn_method = "nndescent"` is specified then many more
  metrics are avaiable, including:

  - `"braycurtis"`

  - `"canberra"`

  - `"chebyshev"`

  - `"dice"`

  - `"hamming"`

  - `"hellinger"`

  - `"jaccard"`

  - `"jensenshannon"`

  - `"kulsinski"`

  - `"rogerstanimoto"`

  - `"russellrao"`

  - `"sokalmichener"`

  - `"sokalsneath"`

  - `"spearmanr"`

  - `"symmetrickl"`

  - `"tsss"`

  - `"yule"`

  For more details see the package documentation of `rnndescent`. For
  `nn_method = "fnn"`, the distance metric is always "euclidean".

  If `X` is a data frame or matrix, then multiple metrics can be
  specified, by passing a list to this argument, where the name of each
  item in the list is one of the metric names above. The value of each
  list item should be a vector giving the names or integer ids of the
  columns to be included in a calculation, e.g.
  `metric = list(euclidean = 1:4, manhattan = 5:10)`.

  Each metric calculation results in a separate fuzzy simplicial set,
  which are intersected together to produce the final set. Metric names
  can be repeated. Because non-numeric columns are removed from the data
  frame, it is safer to use column names than integer ids.

  Factor columns can also be used by specifying the metric name
  `"categorical"`. Factor columns are treated different from numeric
  columns and although multiple factor columns can be specified in a
  vector, each factor column specified is processed individually. If you
  specify a non-factor column, it will be coerced to a factor.

  For a given data block, you may override the `pca` and `pca_center`
  arguments for that block, by providing a list with one unnamed item
  containing the column names or ids, and then any of the `pca` or
  `pca_center` overrides as named items, e.g.
  `metric = list(euclidean = 1:4, manhattan = list(5:10, pca_center = FALSE))`.
  This exists to allow mixed binary and real-valued data to be included
  and to have PCA applied to both, but with centering applied only to
  the real-valued data (it is typical not to apply centering to binary
  data before PCA is applied).

- n_epochs:

  Number of epochs to use during the optimization of the embedded
  coordinates. The default is calculate the number of epochs dynamically
  based on dataset size, to give the same number of edge samples as the
  LargeVis defaults. This is usually substantially larger than the UMAP
  defaults. If `n_epochs = 0`, then coordinates determined by `"init"`
  will be returned.

- learning_rate:

  Initial learning rate used in optimization of the coordinates.

- scale:

  Scaling to apply to `X` if it is a data frame or matrix:

  - `"none"` or `FALSE` or `NULL` No scaling.

  - `"Z"` or `"scale"` or `TRUE` Scale each column to zero mean and
    variance 1.

  - `"maxabs"` Center each column to mean 0, then divide each element by
    the maximum absolute value over the entire matrix.

  - `"range"` Range scale the entire matrix, so the smallest element is
    0 and the largest is 1.

  - `"colrange"` Scale each column in the range (0,1).

  For lvish, the default is `"maxabs"`, for consistency with LargeVis.

- init:

  Type of initialization for the coordinates. Options are:

  - `"spectral"` Spectral embedding using the normalized Laplacian of
    the fuzzy 1-skeleton, with Gaussian noise added.

  - `"normlaplacian"`. Spectral embedding using the normalized Laplacian
    of the fuzzy 1-skeleton, without noise.

  - `"random"`. Coordinates assigned using a uniform random distribution
    between -10 and 10.

  - `"lvrandom"`. Coordinates assigned using a Gaussian distribution
    with standard deviation 1e-4, as used in LargeVis (Tang et
    al., 2016) and t-SNE.

  - `"laplacian"`. Spectral embedding using the Laplacian Eigenmap
    (Belkin and Niyogi, 2002).

  - `"pca"`. The first two principal components from PCA of `X` if `X`
    is a data frame, and from a 2-dimensional classical MDS if `X` is of
    class `"dist"`.

  - `"spca"`. Like `"pca"`, but each dimension is then scaled so the
    standard deviation is 1e-4, to give a distribution similar to that
    used in t-SNE and LargeVis. This is an alias for
    `init = "pca", init_sdev = 1e-4`.

  - `"agspectral"` An "approximate global" modification of `"spectral"`
    which all edges in the graph to a value of 1, and then sets a random
    number of edges (`negative_sample_rate` edges per vertex) to 0.1, to
    approximate the effect of non-local affinities.

  - A matrix of initial coordinates.

  For spectral initializations, (`"spectral"`, `"normlaplacian"`,
  `"laplacian"`, `"agspectral"`), if more than one connected component
  is identified, no spectral initialization is attempted. Instead a
  PCA-based initialization is attempted. If `verbose = TRUE` the number
  of connected components are logged to the console. The existence of
  multiple connected components implies that a global view of the data
  cannot be attained with this initialization. Increasing the value of
  `n_neighbors` may help.

- init_sdev:

  If non-`NULL`, scales each dimension of the initialized coordinates
  (including any user-supplied matrix) to this standard deviation. By
  default no scaling is carried out, except when `init = "spca"`, in
  which case the value is `0.0001`. Scaling the input may help if the
  unscaled versions result in initial coordinates with large inter-point
  distances or outliers. This usually results in small gradients during
  optimization and very little progress being made to the layout.
  Shrinking the initial embedding by rescaling can help under these
  circumstances. Scaling the result of `init = "pca"` is usually
  recommended and `init = "spca"` as an alias for
  `init = "pca", init_sdev = 1e-4` but for the spectral initializations
  the scaled versions usually aren't necessary unless you are using a
  large value of `n_neighbors` (e.g. `n_neighbors = 150` or higher). For
  compatibility with recent versions of the Python UMAP package, if you
  are using `init = "spectral"`, then you should also set
  `init_sdev = "range"`, which will range scale each of the columns
  containing the initial data between 0-10. This is not set by default
  to maintain backwards compatibility with previous versions of uwot.

- repulsion_strength:

  Weighting applied to negative samples in low dimensional embedding
  optimization. Values higher than one will result in greater weight
  being given to negative samples.

- negative_sample_rate:

  The number of negative edge/1-simplex samples to use per positive
  edge/1-simplex sample in optimizing the low dimensional embedding.

- nn_method:

  Method for finding nearest neighbors. Options are:

  - `"fnn"`. Use exact nearest neighbors via the
    [FNN](https://cran.r-project.org/package=FNN) package.

  - `"annoy"` Use approximate nearest neighbors via the
    [RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy) package.

  - `"hnsw"` Use approximate nearest neighbors with the Hierarchical
    Navigable Small World (HNSW) method (Malkov and Yashunin, 2018) via
    the [RcppHNSW](https://cran.r-project.org/package=RcppHNSW) package.
    `RcppHNSW` is not a dependency of this package: this option is only
    available if you have installed `RcppHNSW` yourself. Also, HNSW only
    supports the following arguments for `metric`: `"euclidean"`,
    `"cosine"` and `"correlation"`.

  - `"nndescent"` Use approximate nearest neighbors with the Nearest
    Neighbor Descent method (Dong et al., 2011) via the
    [rnndescent](https://cran.r-project.org/package=rnndescent) package.
    `rnndescent` is not a dependency of this package: this option is
    only available if you have installed `rnndescent` yourself.

  By default, if `X` has less than 4,096 vertices, the exact nearest
  neighbors are found. Otherwise, approximate nearest neighbors are
  used. You may also pass precalculated nearest neighbor data to this
  argument. It must be a list consisting of two elements:

  - `"idx"`. A `n_vertices x n_neighbors` matrix containing the integer
    indexes of the nearest neighbors in `X`. Each vertex is considered
    to be its own nearest neighbor, i.e. `idx[, 1] == 1:n_vertices`.

  - `"dist"`. A `n_vertices x n_neighbors` matrix containing the
    distances of the nearest neighbors.

  Multiple nearest neighbor data (e.g. from two different precomputed
  metrics) can be passed by passing a list containing the nearest
  neighbor data lists as items. The `n_neighbors` parameter is ignored
  when using precomputed nearest neighbor data.

- n_trees:

  Number of trees to build when constructing the nearest neighbor index.
  The more trees specified, the larger the index, but the better the
  results. With `search_k`, determines the accuracy of the Annoy nearest
  neighbor search. Only used if the `nn_method` is `"annoy"`. Sensible
  values are between `10` to `100`.

- search_k:

  Number of nodes to search during the neighbor retrieval. The larger k,
  the more the accurate results, but the longer the search takes. With
  `n_trees`, determines the accuracy of the Annoy nearest neighbor
  search. Only used if the `nn_method` is `"annoy"`.

- n_threads:

  Number of threads to use (except during stochastic gradient descent).
  Default is half the number of concurrent threads supported by the
  system. For nearest neighbor search, only applies if
  `nn_method = "annoy"`. If `n_threads > 1`, then the Annoy index will
  be temporarily written to disk in the location determined by
  [`tempfile`](https://rdrr.io/r/base/tempfile.html).

- n_sgd_threads:

  Number of threads to use during stochastic gradient descent. If set to
  \> 1, then be aware that if `batch = FALSE`, results will *not* be
  reproducible, even if `set.seed` is called with a fixed seed before
  running. Set to `"auto"` to use the same value as `n_threads`.

- grain_size:

  The minimum amount of work to do on each thread. If this value is set
  high enough, then less than `n_threads` or `n_sgd_threads` will be
  used for processing, which might give a performance improvement if the
  overhead of thread management and context switching was outweighing
  the improvement due to concurrent processing. This should be left at
  default (`1`) and work will be spread evenly over all the threads
  specified.

- kernel:

  Type of kernel function to create input probabilities. Can be one of
  `"gauss"` (the default) or `"knn"`. `"gauss"` uses the usual Gaussian
  weighted similarities. `"knn"` assigns equal probabilities to every
  edge in the nearest neighbor graph, and zero otherwise, using
  `perplexity` nearest neighbors. The `n_neighbors` parameter is ignored
  in this case.

- pca:

  If set to a positive integer value, reduce data to this number of
  columns using PCA. Doesn't applied if the distance `metric` is
  `"hamming"`, or the dimensions of the data is larger than the number
  specified (i.e. number of rows and columns must be larger than the
  value of this parameter). If you have \> 100 columns in a data frame
  or matrix, reducing the number of columns in this way may
  substantially increase the performance of the nearest neighbor search
  at the cost of a potential decrease in accuracy. In many t-SNE
  applications, a value of 50 is recommended, although there's no
  guarantee that this is appropriate for all settings.

- pca_center:

  If `TRUE`, center the columns of `X` before carrying out PCA. For
  binary data, it's recommended to set this to `FALSE`.

- pcg_rand:

  If `TRUE`, use the PCG random number generator (O'Neill, 2014) during
  optimization. Otherwise, use the faster (but probably less
  statistically good) Tausworthe "taus88" generator. The default is
  `TRUE`. This parameter has been superseded by `rng_type` – if both are
  set, `rng_type` takes precedence.

- fast_sgd:

  If `TRUE`, then the following combination of parameters is set:
  `pcg_rand = TRUE` and `n_sgd_threads = "auto"`. The default is
  `FALSE`. Setting this to `TRUE` will speed up the stochastic
  optimization phase, but give a potentially less accurate embedding,
  and which will not be exactly reproducible even with a fixed seed. For
  visualization, `fast_sgd = TRUE` will give perfectly good results. For
  more generic dimensionality reduction, it's safer to leave
  `fast_sgd = FALSE`. If `fast_sgd = TRUE`, then user-supplied values of
  `pcg_rand` and `n_sgd_threads`, are ignored.

- ret_nn:

  If `TRUE`, then in addition to the embedding, also return nearest
  neighbor data that can be used as input to `nn_method` to avoid the
  overhead of repeatedly calculating the nearest neighbors when
  manipulating unrelated parameters (e.g. `min_dist`, `n_epochs`,
  `init`). See the "Value" section for the names of the list items. If
  `FALSE`, just return the coordinates. Note that the nearest neighbors
  could be sensitive to data scaling, so be wary of reusing nearest
  neighbor data if modifying the `scale` parameter.

- ret_extra:

  A vector indicating what extra data to return. May contain any
  combination of the following strings:

  - `"nn"` same as setting `ret_nn = TRUE`.

  - `"P"` the high dimensional probability matrix. The graph is returned
    as a sparse symmetric N x N matrix of class
    [dgCMatrix-class](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html),
    where a non-zero entry (i, j) gives the input probability (or
    similarity or affinity) of the edge connecting vertex i and
    vertex j. Note that the graph is further sparsified by removing
    edges with sufficiently low membership strength that they would not
    be sampled by the probabilistic edge sampling employed for
    optimization and therefore the number of non-zero elements in the
    matrix is dependent on `n_epochs`. If you are only interested in the
    fuzzy input graph (e.g. for clustering), setting `n_epochs = 0` will
    avoid any further sparsifying. Be aware that setting
    `binary_edge_weights = TRUE` will affect this graph (all non-zero
    edge weights will be 1).

  - `sigma` a vector of the bandwidths used to calibrate the input
    Gaussians to reproduce the target `"perplexity"`.

- tmpdir:

  Temporary directory to store nearest neighbor indexes during nearest
  neighbor search. Default is
  [`tempdir`](https://rdrr.io/r/base/tempfile.html). The index is only
  written to disk if `n_threads > 1` and `nn_method = "annoy"`;
  otherwise, this parameter is ignored.

- verbose:

  If `TRUE`, log details to the console.

- batch:

  If `TRUE`, then embedding coordinates are updated at the end of each
  epoch rather than during the epoch. In batch mode, results are
  reproducible with a fixed random seed even with `n_sgd_threads > 1`,
  at the cost of a slightly higher memory use. You may also have to
  modify `learning_rate` and increase `n_epochs`, so whether this
  provides a speed increase over the single-threaded optimization is
  likely to be dataset and hardware-dependent.

- opt_args:

  A list of optimizer parameters, used when `batch = TRUE`. The default
  optimization method used is Adam (Kingma and Ba, 2014).

  - `method` The optimization method to use. Either `"adam"` or `"sgd"`
    (stochastic gradient descent). Default: `"adam"`.

  - `beta1` (Adam only). The weighting parameter for the exponential
    moving average of the first moment estimator. Effectively the
    momentum parameter. Should be a floating point value between 0
    and 1. Higher values can smooth oscillatory updates in
    poorly-conditioned situations and may allow for a larger
    `learning_rate` to be specified, but too high can cause divergence.
    Default: `0.5`.

  - `beta2` (Adam only). The weighting parameter for the exponential
    moving average of the uncentered second moment estimator. Should be
    a floating point value between 0 and 1. Controls the degree of
    adaptivity in the step-size. Higher values put more weight on
    previous time steps. Default: `0.9`.

  - `eps` (Adam only). Intended to be a small value to prevent division
    by zero, but in practice can also affect convergence due to its
    interaction with `beta2`. Higher values reduce the effect of the
    step-size adaptivity and bring the behavior closer to stochastic
    gradient descent with momentum. Typical values are between 1e-8 and
    1e-3. Default: `1e-7`.

  - `alpha` The initial learning rate. Default: the value of the
    `learning_rate` parameter.

- epoch_callback:

  A function which will be invoked at the end of every epoch. Its
  signature should be: `(epoch, n_epochs, coords)`, where:

  - `epoch` The current epoch number (between `1` and `n_epochs`).

  - `n_epochs` Number of epochs to use during the optimization of the
    embedded coordinates.

  - `coords` The embedded coordinates as of the end of the current
    epoch, as a matrix with dimensions (N, `n_components`).

- pca_method:

  Method to carry out any PCA dimensionality reduction when the `pca`
  parameter is specified. Allowed values are:

  - `"irlba"`. Uses
    [`prcomp_irlba`](https://rdrr.io/pkg/irlba/man/prcomp_irlba.html)
    from the [irlba](https://cran.r-project.org/package=irlba) package.

  - `"rsvd"`. Uses 5 iterations of
    [`svdr`](https://rdrr.io/pkg/irlba/man/svdr.html) from the
    [irlba](https://cran.r-project.org/package=irlba) package. This is
    likely to give much faster but potentially less accurate results
    than using `"irlba"`. For the purposes of nearest neighbor
    calculation and coordinates initialization, any loss of accuracy
    doesn't seem to matter much.

  - `"bigstatsr"`. Uses
    [`big_randomSVD`](https://privefl.github.io/bigstatsr/reference/big_randomSVD.html)
    from the [bigstatsr](https://cran.r-project.org/package=bigstatsr)
    package. The SVD methods used in `bigstatsr` may be faster on
    systems without access to efficient linear algebra libraries (e.g.
    Windows). **Note**: `bigstatsr` is *not* a dependency of uwot: if
    you choose to use this package for PCA, you *must* install it
    yourself.

  - `"svd"`. Uses [`svd`](https://rdrr.io/r/base/svd.html) for the SVD.
    This is likely to be slow for all but the smallest datasets.

  - `"auto"` (the default). Uses `"irlba"`, unless more than 50 case
    `"svd"` is used.

- binary_edge_weights:

  If `TRUE` then edge weights in the input graph are treated as binary
  (0/1) rather than real valued. This affects the sampling frequency of
  neighbors and is the strategy used by the PaCMAP method (Wang and
  co-workers, 2020). Practical (Böhm and co-workers, 2020) and
  theoretical (Damrich and Hamprecht, 2021) work suggests this has
  little effect on UMAP's performance.

- nn_args:

  A list containing additional arguments to pass to the nearest neighbor
  method. For `nn_method = "annoy"`, you can specify `"n_trees"` and
  `"search_k"`, and these will override the `n_trees` and `search_k`
  parameters. For `nn_method = "hnsw"`, you may specify the following
  arguments:

  - `M` The maximum number of neighbors to keep for each vertex.
    Reasonable values are `2` to `100`. Higher values give better recall
    at the cost of more memory. Default value is `16`.

  - `ef_construction` A positive integer specifying the size of the
    dynamic list used during index construction. A higher value will
    provide better results at the cost of a longer time to build the
    index. Default is `200`.

  - `ef` A positive integer specifying the size of the dynamic list used
    during search. This cannot be smaller than `n_neighbors` and cannot
    be higher than the number of items in the index. Default is `10`.

  For `nn_method = "nndescent"`, you may specify the following
  arguments:

  - `n_trees` The number of trees to use in a random projection forest
    to initialize the search. A larger number will give more accurate
    results at the cost of a longer computation time. The default of
    `NULL` means that the number is chosen based on the number of
    observations in `X`.

  - `max_candidates` The number of potential neighbors to explore per
    iteration. By default, this is set to `n_neighbors` or `60`,
    whichever is smaller. A larger number will give more accurate
    results at the cost of a longer computation time.

  - `n_iters` The number of iterations to run the search. A larger
    number will give more accurate results at the cost of a longer
    computation time. By default, this will be chosen based on the
    number of observations in `X`. You may also need to modify the
    convergence criterion `delta`.

  - `delta` The minimum relative change in the neighbor graph allowed
    before early stopping. Should be a value between 0 and 1. The
    smaller the value, the smaller the amount of progress between
    iterations is allowed. Default value of `0.001` means that at least
    0.1 neighbor graph must be updated at each iteration.

  - `init` How to initialize the nearest neighbor descent. By default
    this is set to `"tree"` and uses a random project forest. If you set
    this to `"rand"`, then a random selection is used. Usually this is
    less accurate than using RP trees, but for high-dimensional cases,
    there may be little difference in the quality of the initialization
    and random initialization will be a lot faster. If you set this to
    `"rand"`, then the `n_trees` parameter is ignored.

- rng_type:

  The type of random number generator to use during optimization. One
  of:

  - `"pcg"`. Use the PCG random number generator (O'Neill, 2014).

  - `"tausworthe"`. Use the Tausworthe "taus88" generator.

  - `"deterministic"`. Use a deterministic number generator. This isn't
    actually random, but may provide enough variation in the negative
    sampling to give a good embedding and can provide a noticeable
    speed-up.

  For backwards compatibility, by default this is unset and the choice
  of `pcg_rand` is used (making "pcg" the effective default).

## Value

A matrix of optimized coordinates, or:

- if `ret_nn = TRUE` (or `ret_extra` contains `"nn"`), returns the
  nearest neighbor data as a list called `nn`. This contains one list
  for each `metric` calculated, itself containing a matrix `idx` with
  the integer ids of the neighbors; and a matrix `dist` with the
  distances. The `nn` list (or a sub-list) can be used as input to the
  `nn_method` parameter.

- if `ret_extra` contains `"P"`, returns the high dimensional
  probability matrix as a sparse matrix called `P`, of type
  [dgCMatrix-class](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html).

- if `ret_extra` contains `"sigma"`, returns a vector of the high
  dimensional gaussian bandwidths for each point, and `"dint"` a vector
  of estimates of the intrinsic dimensionality at each point, based on
  the method given by Lee and co-workers (2015).

The returned list contains the combined data from any combination of
specifying `ret_nn` and `ret_extra`.

## Details

`lvish` differs from the official LargeVis implementation in the
following:

- Only the nearest-neighbor index search phase is multi-threaded.

- Matrix input data is not normalized.

- The `n_trees` parameter cannot be dynamically chosen based on data set
  size.

- Nearest neighbor results are not refined via the
  neighbor-of-my-neighbor method. The `search_k` parameter is twice as
  large than default to compensate.

- Gradient values are clipped to `4.0` rather than `5.0`.

- Negative edges are generated by uniform sampling of vertexes rather
  than their degree ^ 0.75.

- The default number of samples is much reduced. The default number of
  epochs, `n_epochs`, is set to `5000`, much larger than for
  [`umap`](https://jlmelville.github.io/uwot/reference/umap.md), but may
  need to be increased further depending on your dataset. Using
  `init = "spectral"` can help.

## References

Belkin, M., & Niyogi, P. (2002). Laplacian eigenmaps and spectral
techniques for embedding and clustering. In *Advances in neural
information processing systems* (pp. 585-591).
<http://papers.nips.cc/paper/1961-laplacian-eigenmaps-and-spectral-techniques-for-embedding-and-clustering.pdf>

Böhm, J. N., Berens, P., & Kobak, D. (2020). A unifying perspective on
neighbor embeddings along the attraction-repulsion spectrum. *arXiv
preprint* *arXiv:2007.08902*. <https://arxiv.org/abs/2007.08902>

Damrich, S., & Hamprecht, F. A. (2021). On UMAP's true loss function.
*Advances in Neural Information Processing Systems*, *34*.
<https://proceedings.neurips.cc/paper/2021/hash/2de5d16682c3c35007e4e92982f1a2ba-Abstract.html>

Dong, W., Moses, C., & Li, K. (2011, March). Efficient k-nearest
neighbor graph construction for generic similarity measures. In
*Proceedings of the 20th international conference on World Wide Web*
(pp. 577-586). ACM.
[doi:10.1145/1963405.1963487](https://doi.org/10.1145/1963405.1963487) .

Kingma, D. P., & Ba, J. (2014). Adam: A method for stochastic
optimization. *arXiv preprint* *arXiv*:1412.6980.
<https://arxiv.org/abs/1412.6980>

Lee, J. A., Peluffo-Ordóñez, D. H., & Verleysen, M. (2015). Multi-scale
similarities in stochastic neighbour embedding: Reducing dimensionality
while preserving both local and global structure. *Neurocomputing*,
*169*, 246-261.

Malkov, Y. A., & Yashunin, D. A. (2018). Efficient and robust
approximate nearest neighbor search using hierarchical navigable small
world graphs. *IEEE transactions on pattern analysis and machine
intelligence*, *42*(4), 824-836.

McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform Manifold
Approximation and Projection for Dimension Reduction *arXiv preprint*
*arXiv*:1802.03426. <https://arxiv.org/abs/1802.03426>

O'Neill, M. E. (2014). *PCG: A family of simple fast space-efficient
statistically good algorithms for random number generation* (Report No.
HMC-CS-2014-0905). Harvey Mudd College.

Tang, J., Liu, J., Zhang, M., & Mei, Q. (2016, April). Visualizing
large-scale and high-dimensional data. In *Proceedings of the 25th
International Conference on World Wide Web* (pp. 287-297). International
World Wide Web Conferences Steering Committee.
<https://arxiv.org/abs/1602.00370>

Van der Maaten, L., & Hinton, G. (2008). Visualizing data using t-SNE.
*Journal of Machine Learning Research*, *9* (2579-2605).
<https://www.jmlr.org/papers/v9/vandermaaten08a.html>

Wang, Y., Huang, H., Rudin, C., & Shaposhnik, Y. (2021). Understanding
How Dimension Reduction Tools Work: An Empirical Approach to Deciphering
t-SNE, UMAP, TriMap, and PaCMAP for Data Visualization. *Journal of
Machine Learning Research*, *22*(201), 1-73.
<https://www.jmlr.org/papers/v22/20-1061.html>

## Examples

``` r
# Default number of epochs is much larger than for UMAP, assumes random
# initialization. Use perplexity rather than n_neighbors to control the size
# of the local neighborhood 20 epochs may be too small for a random
# initialization
iris_lvish <- lvish(iris,
  perplexity = 50, learning_rate = 0.5,
  init = "random", n_epochs = 20
)
```
