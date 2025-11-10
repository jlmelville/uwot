# Dimensionality Reduction Using t-Distributed UMAP (t-UMAP)

A faster (but less flexible) version of the UMAP (McInnes et al, 2018)
gradient. For more detail on UMAP, see the
[`umap`](https://jlmelville.github.io/uwot/reference/umap.md) function.

## Usage

``` r
tumap(
  X,
  n_neighbors = 15,
  n_components = 2,
  metric = "euclidean",
  n_epochs = NULL,
  learning_rate = 1,
  scale = FALSE,
  init = "spectral",
  init_sdev = NULL,
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  bandwidth = 1,
  repulsion_strength = 1,
  negative_sample_rate = 5,
  nn_method = NULL,
  n_trees = 50,
  search_k = 2 * n_neighbors * n_trees,
  n_threads = NULL,
  n_sgd_threads = 0,
  grain_size = 1,
  y = NULL,
  target_n_neighbors = n_neighbors,
  target_metric = "euclidean",
  target_weight = 0.5,
  pca = NULL,
  pca_center = TRUE,
  pcg_rand = TRUE,
  fast_sgd = FALSE,
  ret_model = FALSE,
  ret_nn = FALSE,
  ret_extra = c(),
  tmpdir = tempdir(),
  verbose = getOption("verbose", TRUE),
  batch = FALSE,
  opt_args = NULL,
  epoch_callback = NULL,
  pca_method = NULL,
  binary_edge_weights = FALSE,
  seed = NULL,
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

- n_neighbors:

  The size of local neighborhood (in terms of number of neighboring
  sample points) used for manifold approximation. Larger values result
  in more global views of the manifold, while smaller values result in
  more local data being preserved. In general values should be in the
  range `2` to `100`.

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
  coordinates. By default, this value is set to `500` for datasets
  containing 10,000 vertices or less, and `200` otherwise. If
  `n_epochs = 0`, then coordinates determined by `"init"` will be
  returned.

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

  For t-UMAP, the default is `"none"`.

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
    used in t-SNE. This is an alias for
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

- set_op_mix_ratio:

  Interpolate between (fuzzy) union and intersection as the set
  operation used to combine local fuzzy simplicial sets to obtain a
  global fuzzy simplicial sets. Both fuzzy set operations use the
  product t-norm. The value of this parameter should be between `0.0`
  and `1.0`; a value of `1.0` will use a pure fuzzy union, while `0.0`
  will use a pure fuzzy intersection.

- local_connectivity:

  The local connectivity required – i.e. the number of nearest neighbors
  that should be assumed to be connected at a local level. The higher
  this value the more connected the manifold becomes locally. In
  practice this should be not more than the local intrinsic dimension of
  the manifold.

- bandwidth:

  The effective bandwidth of the kernel if we view the algorithm as
  similar to Laplacian Eigenmaps. Larger values induce more connectivity
  and a more global view of the data, smaller values concentrate more
  locally.

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
    supports the following arguments for `metric` and `target_metric`:
    `"euclidean"`, `"cosine"` and `"correlation"`.

  - `"nndescent"` Use approximate nearest neighbors with the Nearest
    Neighbor Descent method (Dong et al., 2011) via the
    [rnndescent](https://cran.r-project.org/package=rnndescent) package.
    `rnndescent` is not a dependency of this package: this option is
    only available if you have installed `rnndescent` yourself.

  By default, if `X` has less than 4,096 vertices, the exact nearest
  neighbors are found. Otherwise, approximate nearest neighbors are
  used. You may also pass pre-calculated nearest neighbor data to this
  argument. It must be one of two formats, either a list consisting of
  two elements:

  - `"idx"`. A `n_vertices x n_neighbors` matrix containing the integer
    indexes of the nearest neighbors in `X`. Each vertex is considered
    to be its own nearest neighbor, i.e. `idx[, 1] == 1:n_vertices`.

  - `"dist"`. A `n_vertices x n_neighbors` matrix containing the
    distances of the nearest neighbors.

  or a sparse distance matrix of type `dgCMatrix`, with dimensions
  `n_vertices x n_vertices`. Distances should be arranged by column,
  i.e. a non-zero entry in row `j` of the `i`th column indicates that
  the `j`th observation in `X` is a nearest neighbor of the `i`th
  observation with the distance given by the value of that element. The
  `n_neighbors` parameter is ignored when using precomputed nearest
  neighbor data. If using the sparse distance matrix input, each column
  can contain a different number of neighbors.

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

- y:

  Optional target data for supervised dimension reduction. Can be a
  vector, matrix or data frame. Use the `target_metric` parameter to
  specify the metrics to use, using the same syntax as `metric`. Usually
  either a single numeric or factor column is used, but more complex
  formats are possible. The following types are allowed:

  - Factor columns with the same length as `X`. `NA` is allowed for any
    observation with an unknown level, in which case UMAP operates as a
    form of semi-supervised learning. Each column is treated separately.

  - Numeric data. `NA` is *not* allowed in this case. Use the parameter
    `target_n_neighbors` to set the number of neighbors used with `y`.
    If unset, `n_neighbors` is used. Unlike factors, numeric columns are
    grouped into one block unless `target_metric` specifies otherwise.
    For example, if you wish columns `a` and `b` to be treated
    separately, specify
    `target_metric = list(euclidean = "a", euclidean = "b")`. Otherwise,
    the data will be effectively treated as a matrix with two columns.

  - Nearest neighbor data, consisting of a list of two matrices, `idx`
    and `dist`. These represent the precalculated nearest neighbor
    indices and distances, respectively. This is the same format as that
    expected for precalculated data in `nn_method`. This format assumes
    that the underlying data was a numeric vector. Any user-supplied
    value of the `target_n_neighbors` parameter is ignored in this case,
    because the the number of columns in the matrices is used for the
    value. Multiple nearest neighbor data using different metrics can be
    supplied by passing a list of these lists.

  Unlike `X`, all factor columns included in `y` are automatically used.

- target_n_neighbors:

  Number of nearest neighbors to use to construct the target simplicial
  set. Default value is `n_neighbors`. Applies only if `y` is non-`NULL`
  and `numeric`.

- target_metric:

  The metric used to measure distance for `y` if using supervised
  dimension reduction. Used only if `y` is numeric.

- target_weight:

  Weighting factor between data topology and target topology. A value of
  0.0 weights entirely on data, a value of 1.0 weights entirely on
  target. The default of 0.5 balances the weighting equally between data
  and target. Only applies if `y` is non-`NULL`.

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

- ret_model:

  If `TRUE`, then return extra data that can be used to add new data to
  an existing embedding via
  [`umap_transform`](https://jlmelville.github.io/uwot/reference/umap_transform.md).
  The embedded coordinates are returned as the list item `embedding`. If
  `FALSE`, just return the coordinates. This parameter can be used in
  conjunction with `ret_nn` and `ret_extra`. Note that some settings are
  incompatible with the production of a UMAP model: external neighbor
  data (passed via a list to `nn_method`), and factor columns that were
  included via the `metric` parameter. In the latter case, the model
  produced is based only on the numeric data. A transformation using new
  data is possible, but the factor columns in the new data are ignored.
  Note that setting `ret_model = TRUE` forces the use of the approximate
  nearest neighbors method. Because small datasets would otherwise use
  exact nearest neighbor calculations, setting `ret_model = TRUE` means
  that different results may be returned for small datasets in terms of
  both the returned nearest neighbors (if requested) and the final
  embedded coordinates, compared to `ret_model = FALSE`, even if the
  random number seed is fixed. To avoid this, explicitly set
  `nn_method = "annoy"` in the `ret_model = FALSE` case.

- ret_nn:

  If `TRUE`, then in addition to the embedding, also return nearest
  neighbor data that can be used as input to `nn_method` to avoid the
  overhead of repeatedly calculating the nearest neighbors when
  manipulating unrelated parameters (e.g. `min_dist`, `n_epochs`,
  `init`). See the "Value" section for the names of the list items. If
  `FALSE`, just return the coordinates. Note that the nearest neighbors
  could be sensitive to data scaling, so be wary of reusing nearest
  neighbor data if modifying the `scale` parameter. This parameter can
  be used in conjunction with `ret_model` and `ret_extra`.

- ret_extra:

  A vector indicating what extra data to return. May contain any
  combination of the following strings:

  - `"model"` Same as setting `ret_model = TRUE`.

  - `"nn"` Same as setting `ret_nn = TRUE`.

  - `"fgraph"` the high dimensional fuzzy graph (i.e. the fuzzy
    simplicial set of the merged local views of the input data). The
    graph is returned as a sparse symmetric N x N matrix of class
    [dgCMatrix-class](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html),
    where a non-zero entry (i, j) gives the membership strength of the
    edge connecting vertex i and vertex j. This can be considered
    analogous to the input probability (or similarity or affinity) used
    in t-SNE and LargeVis. Note that the graph is further sparsified by
    removing edges with sufficiently low membership strength that they
    would not be sampled by the probabilistic edge sampling employed for
    optimization and therefore the number of non-zero elements in the
    matrix is dependent on `n_epochs`. If you are only interested in the
    fuzzy input graph (e.g. for clustering), setting `n_epochs = 0` will
    avoid any further sparsifying. Be aware that setting
    `binary_edge_weights = TRUE` will affect this graph (all non-zero
    edge weights will be 1).

  - `"sigma"` the normalization value for each observation in the
    dataset when constructing the smoothed distances to each of its
    neighbors. This gives some sense of the local density of each
    observation in the high dimensional space: higher values of `sigma`
    indicate a higher dispersion or lower density.

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

- seed:

  Integer seed to use to initialize the random number generator state.
  Combined with `n_sgd_threads = 1` or `batch = TRUE`, this should give
  consistent output across multiple runs on a given installation.
  Setting this value is equivalent to calling
  [`set.seed`](https://rdrr.io/r/base/Random.html), but it may be more
  convenient in some situations than having to call a separate function.
  The default is to not set a seed. If `ret_model = TRUE`, the seed will
  be stored in the output model and then used to set the seed inside
  [`umap_transform`](https://jlmelville.github.io/uwot/reference/umap_transform.md).

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

  - `pruning_degree_multiplier` The maximum number of edges per node to
    retain in the search graph, relative to `n_neighbors`. A larger
    value will give more accurate results at the cost of a longer
    computation time. Default is `1.5`. This parameter only affects
    neighbor search when transforming new data with
    [`umap_transform`](https://jlmelville.github.io/uwot/reference/umap_transform.md).

  - `epsilon` Controls the degree of the back-tracking when traversing
    the search graph. Setting this to `0.0` will do a greedy search with
    no back-tracking. A larger value will give more accurate results at
    the cost of a longer computation time. Default is `0.1`. This
    parameter only affects neighbor search when transforming new data
    with
    [`umap_transform`](https://jlmelville.github.io/uwot/reference/umap_transform.md).

  - `max_search_fraction` Specifies the maximum fraction of the search
    graph to traverse. By default, this is set to `1.0`, so the entire
    graph (i.e. all items in `X`) may be visited. You may want to set
    this to a smaller value if you have a very large dataset (in
    conjunction with `epsilon`) to avoid an inefficient exhaustive
    search of the data in `X`. This parameter only affects neighbor
    search when transforming new data with
    [`umap_transform`](https://jlmelville.github.io/uwot/reference/umap_transform.md).

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

  - `pruning_degree_multiplier` The maximum number of edges per node to
    retain in the search graph, relative to `n_neighbors`. A larger
    value will give more accurate results at the cost of a longer
    computation time. Default is `1.5`. This parameter only affects
    neighbor search when transforming new data with
    [`umap_transform`](https://jlmelville.github.io/uwot/reference/umap_transform.md).

  - `epsilon` Controls the degree of the back-tracking when traversing
    the search graph. Setting this to `0.0` will do a greedy search with
    no back-tracking. A larger value will give more accurate results at
    the cost of a longer computation time. Default is `0.1`. This
    parameter only affects neighbor search when transforming new data
    with
    [`umap_transform`](https://jlmelville.github.io/uwot/reference/umap_transform.md).

  - `max_search_fraction` Specifies the maximum fraction of the search
    graph to traverse. By default, this is set to `1.0`, so the entire
    graph (i.e. all items in `X`) may be visited. You may want to set
    this to a smaller value if you have a very large dataset (in
    conjunction with `epsilon`) to avoid an inefficient exhaustive
    search of the data in `X`. This parameter only affects neighbor
    search when transforming new data with
    [`umap_transform`](https://jlmelville.github.io/uwot/reference/umap_transform.md).

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

- if `ret_model = TRUE` (or `ret_extra` contains `"model"`), returns a
  list containing extra information that can be used to add new data to
  an existing embedding via
  [`umap_transform`](https://jlmelville.github.io/uwot/reference/umap_transform.md).
  In this case, the coordinates are available in the list item
  `embedding`. **NOTE**: The contents of the `model` list should *not*
  be considered stable or part of the public API, and are purposely left
  undocumented.

- if `ret_nn = TRUE` (or `ret_extra` contains `"nn"`), returns the
  nearest neighbor data as a list called `nn`. This contains one list
  for each `metric` calculated, itself containing a matrix `idx` with
  the integer ids of the neighbors; and a matrix `dist` with the
  distances. The `nn` list (or a sub-list) can be used as input to the
  `nn_method` parameter.

- if `ret_extra` contains `"fgraph"` returns the high dimensional fuzzy
  graph as a sparse matrix called `fgraph`, of type
  [dgCMatrix-class](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html).

- if `ret_extra` contains `"sigma"`, returns a vector of the smooth knn
  distance normalization terms for each observation as `"sigma"` and a
  vector `"rho"` containing the largest distance to the locally
  connected neighbors of each observation.

- if `ret_extra` contains `"localr"`, returns a vector of the estimated
  local radii, the sum of `"sigma"` and `"rho"`.

The returned list contains the combined data from any combination of
specifying `ret_model`, `ret_nn` and `ret_extra`.

## Details

By setting the UMAP curve parameters `a` and `b` to `1`, you get back
the Cauchy distribution as used in t-SNE (van der Maaten and Hinton,
2008) and LargeVis (Tang et al., 2016). It also results in a
substantially simplified gradient expression. This can give a speed
improvement of around 50%.

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
iris_tumap <- tumap(iris, n_neighbors = 50, learning_rate = 0.5)
```
