# Similarity Graph

Create a graph (as a sparse symmetric weighted adjacency matrix)
representing the similarities between items in a data set. No
dimensionality reduction is carried out. By default, the similarities
are calculated using the merged fuzzy simplicial set approach in the
Uniform Manifold Approximation and Projection (UMAP) method (McInnes et
al., 2018), but the approach from LargeVis (Tang et al., 2016) can also
be used.

## Usage

``` r
similarity_graph(
  X = NULL,
  n_neighbors = NULL,
  metric = "euclidean",
  scale = NULL,
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  nn_method = NULL,
  n_trees = 50,
  search_k = 2 * n_neighbors * n_trees,
  perplexity = 50,
  method = "umap",
  y = NULL,
  target_n_neighbors = n_neighbors,
  target_metric = "euclidean",
  target_weight = 0.5,
  pca = NULL,
  pca_center = TRUE,
  ret_extra = c(),
  n_threads = NULL,
  grain_size = 1,
  kernel = "gauss",
  tmpdir = tempdir(),
  verbose = getOption("verbose", TRUE),
  pca_method = NULL,
  binary_edge_weights = FALSE,
  nn_args = list()
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
  neighbor data is passed to `nn_method`.

- n_neighbors:

  The size of local neighborhood (in terms of number of neighboring
  sample points) used for manifold approximation. Larger values result
  in more global views of the manifold, while smaller values result in
  more local data being preserved. In general values should be in the
  range `2` to `100`.

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

  For `method` `"umap"`, the default is `"none"`. For `"largevis"`, the
  default is `"maxabs"`.

- set_op_mix_ratio:

  Interpolate between (fuzzy) union and intersection as the set
  operation used to combine local fuzzy simplicial sets to obtain a
  global fuzzy simplicial sets. Both fuzzy set operations use the
  product t-norm. The value of this parameter should be between `0.0`
  and `1.0`; a value of `1.0` will use a pure fuzzy union, while `0.0`
  will use a pure fuzzy intersection. Ignored if `method = "largevis"`

- local_connectivity:

  The local connectivity required – i.e. the number of nearest neighbors
  that should be assumed to be connected at a local level. The higher
  this value the more connected the manifold becomes locally. In
  practice this should be not more than the local intrinsic dimension of
  the manifold. Ignored if `method = "largevis"`.

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

- perplexity:

  Used only if `method = "largevis"`. Controls the size of the local
  neighborhood used for manifold approximation. Should be a value
  between 1 and one less than the number of items in `X`. If specified,
  you should *not* specify a value for `n_neighbors` unless you know
  what you are doing.

- method:

  How to generate the similarities between items. One of:

  - `"umap"` The UMAP method of McInnes et al. (2018).

  - `"largevis"` The LargeVis method of Tang et al. (2016).

- y:

  Optional target data to add supervised or semi-supervised weighting to
  the similarity graph . Can be a vector, matrix or data frame. Use the
  `target_metric` parameter to specify the metrics to use, using the
  same syntax as `metric`. Usually either a single numeric or factor
  column is used, but more complex formats are possible. The following
  types are allowed:

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
  This parameter is ignored if `method = "largevis"`.

- target_n_neighbors:

  Number of nearest neighbors to use to construct the target simplicial
  set. Default value is `n_neighbors`. Applies only if `y` is non-`NULL`
  and `numeric`. This parameter is ignored if `method = "largevis"`.

- target_metric:

  The metric used to measure distance for `y` if using supervised
  dimension reduction. Used only if `y` is numeric. This parameter is
  ignored if `method = "largevis"`.

- target_weight:

  Weighting factor between data topology and target topology. A value of
  0.0 weights entirely on data, a value of 1.0 weights entirely on
  target. The default of 0.5 balances the weighting equally between data
  and target. Only applies if `y` is non-`NULL`. This parameter is
  ignored if `method = "largevis"`.

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

- ret_extra:

  A vector indicating what extra data to return. May contain any
  combination of the following strings:

  - `"nn"` nearest neighbor data that can be used as input to
    `nn_method` to avoid the overhead of repeatedly calculating the
    nearest neighbors when manipulating unrelated parameters. See the
    "Value" section for the names of the list items. Note that the
    nearest neighbors could be sensitive to data scaling, so be wary of
    reusing nearest neighbor data if modifying the `scale` parameter.

  - `"sigma"` the normalization value for each observation in the
    dataset when constructing the smoothed distances to each of its
    neighbors. This gives some sense of the local density of each
    observation in the high dimensional space: higher values of `sigma`
    indicate a higher dispersion or lower density.

- n_threads:

  Number of threads to use. Default is half the number of concurrent
  threads supported by the system. For nearest neighbor search, only
  applies if `nn_method = "annoy"`. If `n_threads > 1`, then the Annoy
  index will be temporarily written to disk in the location determined
  by [`tempfile`](https://rdrr.io/r/base/tempfile.html).

- grain_size:

  The minimum amount of work to do on each thread. If this value is set
  high enough, then less than `n_threads` will be used for processing,
  which might give a performance improvement if the overhead of thread
  management and context switching was outweighing the improvement due
  to concurrent processing. This should be left at default (`1`) and
  work will be spread evenly over all the threads specified.

- kernel:

  Used only if `method = "largevis"`. Type of kernel function to create
  input similiarties. Can be one of `"gauss"` (the default) or `"knn"`.
  `"gauss"` uses the usual Gaussian weighted similarities. `"knn"`
  assigns equal similiarties. to every edge in the nearest neighbor
  graph, and zero otherwise, using `perplexity` nearest neighbors. The
  `n_neighbors` parameter is ignored in this case.

- tmpdir:

  Temporary directory to store nearest neighbor indexes during nearest
  neighbor search. Default is
  [`tempdir`](https://rdrr.io/r/base/tempfile.html). The index is only
  written to disk if `n_threads > 1` and `nn_method = "annoy"`;
  otherwise, this parameter is ignored.

- verbose:

  If `TRUE`, log details to the console.

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

  If `TRUE` then edge weights of the returned graph are binary (0/1)
  rather than reflecting the degree of similarity.

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

## Value

A sparse symmetrized matrix of the similarities between the items in `X`
or if `nn_method` contains pre-computed nearest neighbor data, the items
in `nn_method`. Because of the symmetrization, there may be more
non-zero items in each column than the specified value of `n_neighbors`
(or pre-computed neighbors in `nn_method`). If `ret_extra` is specified
then the return value will be a list containing:

- `similarity_graph` the similarity graph as a sparse matrix as
  described above.

- `nn` (if `ret_extra` contained `"nn"`) the nearest neighbor data as a
  list called `nn`. This contains one list for each `metric` calculated,
  itself containing a matrix `idx` with the integer ids of the
  neighbors; and a matrix `dist` with the distances. The `nn` list (or a
  sub-list) can be used as input to the `nn_method` parameter.

- `sigma` (if `ret_extra` contains `"sigma"`), a vector of calibrated
  parameters, one for each item in the input data, reflecting the local
  data density for that item. The exact definition of the values depends
  on the choice of the `method` parameter.

- `rho` (if `ret_extra` contains `"sigma"`), a vector containing the
  largest distance to the locally connected neighbors of each item in
  the input data. This will exist only if `method = "umap"`.

- `localr` (if `ret_extra` contains `"localr"`) a vector of the
  estimated local radii, the sum of `"sigma"` and `"rho"`. This will
  exist only if `method = "umap"`.

## Details

This is equivalent to running
[`umap`](https://jlmelville.github.io/uwot/reference/umap.md) with the
`ret_extra = c("fgraph")` parameter, but without the overhead of
calculating (or returning) the optimized low-dimensional coordinates.

## References

Dong, W., Moses, C., & Li, K. (2011, March). Efficient k-nearest
neighbor graph construction for generic similarity measures. In
*Proceedings of the 20th international conference on World Wide Web*
(pp. 577-586). ACM.
[doi:10.1145/1963405.1963487](https://doi.org/10.1145/1963405.1963487) .

Malkov, Y. A., & Yashunin, D. A. (2018). Efficient and robust
approximate nearest neighbor search using hierarchical navigable small
world graphs. *IEEE transactions on pattern analysis and machine
intelligence*, *42*(4), 824-836.

McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform Manifold
Approximation and Projection for Dimension Reduction *arXiv preprint*
*arXiv*:1802.03426. <https://arxiv.org/abs/1802.03426>

Tang, J., Liu, J., Zhang, M., & Mei, Q. (2016, April). Visualizing
large-scale and high-dimensional data. In *Proceedings of the 25th
International Conference on World Wide Web* (pp. 287-297). International
World Wide Web Conferences Steering Committee.
<https://arxiv.org/abs/1602.00370>

## Examples

``` r
iris30 <- iris[c(1:10, 51:60, 101:110), ]

# return a 30 x 30 sparse matrix with similarity data based on 10 nearest
# neighbors per item
iris30_sim_graph <- similarity_graph(iris30, n_neighbors = 10)

# Default is to use the UMAP method of calculating similarities, but LargeVis
# is also available: for that method, use perplexity instead of n_neighbors
# to control neighborhood size. Use ret_extra = "nn" to return nearest
# neighbor data as well as the similarity graph. Return value is a list
# containing similarity_graph' and 'nn' items.
iris30_lv_graph <- similarity_graph(iris30,
  perplexity = 10,
  method = "largevis", ret_extra = "nn"
)
# If you have the neighbor information you don't need the original data
iris30_lv_graph_nn <- similarity_graph(
  nn_method = iris30_lv_graph$nn,
  perplexity = 10, method = "largevis"
)
all(iris30_lv_graph_nn == iris30_lv_graph$similarity_graph)
#> [1] TRUE
```
