# Optimize Graph Layout

Carry out dimensionality reduction on an input graph, where the
distances in the low dimensional space attempt to reproduce the neighbor
relations in the input data. By default, the cost function used to
optimize the output coordinates use the Uniform Manifold Approximation
and Projection (UMAP) method (McInnes et al., 2018), but the approach
from LargeVis (Tang et al., 2016) can also be used. This function can be
used to produce a low dimensional representation of the graph produced
by
[`similarity_graph`](https://jlmelville.github.io/uwot/reference/similarity_graph.md).

## Usage

``` r
optimize_graph_layout(
  graph,
  X = NULL,
  n_components = 2,
  n_epochs = NULL,
  learning_rate = 1,
  init = "spectral",
  init_sdev = NULL,
  spread = 1,
  min_dist = 0.01,
  repulsion_strength = 1,
  negative_sample_rate = 5,
  a = NULL,
  b = NULL,
  method = "umap",
  approx_pow = FALSE,
  pcg_rand = TRUE,
  fast_sgd = FALSE,
  n_sgd_threads = 0,
  grain_size = 1,
  verbose = getOption("verbose", TRUE),
  batch = FALSE,
  opt_args = NULL,
  epoch_callback = NULL,
  pca_method = NULL,
  binary_edge_weights = FALSE,
  rng_type = NULL
)
```

## Arguments

- graph:

  A sparse, symmetric N x N weighted adjacency matrix representing a
  graph. Non-zero entries indicate an edge between two nodes with a
  given edge weight. There can be a varying number of non-zero entries
  in each row/column.

- X:

  Optional input data. Used only for PCA-based initialization.

- n_components:

  The dimension of the space to embed into. This defaults to `2` to
  provide easy visualization, but can reasonably be set to any integer
  value in the range `2` to `100`.

- n_epochs:

  Number of epochs to use during the optimization of the embedded
  coordinates. By default, this value is set to `500` for datasets
  containing 10,000 vertices or less, and `200` otherwise. If
  `n_epochs = 0`, then coordinates determined by `"init"` will be
  returned. For UMAP, the default is `"none"`.

- learning_rate:

  Initial learning rate used in optimization of the coordinates.

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

  - `"laplacian"`. Spectral embedding using the Laplacian Eigenmap.

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

- spread:

  The effective scale of embedded points. In combination with
  `min_dist`, this determines how clustered/clumped the embedded points
  are.

- min_dist:

  The effective minimum distance between embedded points. Smaller values
  will result in a more clustered/clumped embedding where nearby points
  on the manifold are drawn closer together, while larger values will
  result on a more even dispersal of points. The value should be set
  relative to the `spread` value, which determines the scale at which
  embedded points will be spread out.

- repulsion_strength:

  Weighting applied to negative samples in low dimensional embedding
  optimization. Values higher than one will result in greater weight
  being given to negative samples.

- negative_sample_rate:

  The number of negative edge/1-simplex samples to use per positive
  edge/1-simplex sample in optimizing the low dimensional embedding.

- a:

  More specific parameters controlling the embedding. If `NULL` these
  values are set automatically as determined by `min_dist` and `spread`.

- b:

  More specific parameters controlling the embedding. If `NULL` these
  values are set automatically as determined by `min_dist` and `spread`.

- method:

  Cost function to optimize. One of:

  - `"umap"`. The UMAP method of McInnes and co-workers (2018).

  - `"tumap"`. UMAP with the `a` and `b` parameters fixed to 1.

  - `"largevis"`. The LargeVis method Tang and co-workers (2016).

- approx_pow:

  If `TRUE`, use an approximation to the power function in the UMAP
  gradient, from
  <https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/>.

- pcg_rand:

  If `TRUE`, use the PCG random number generator (O'Neill, 2014) during
  optimization. Otherwise, use the faster (but probably less
  statistically good) Tausworthe "taus88" generator. The default is
  `TRUE`. This parameter has been superseded by `rng_type` – if both are
  set, `rng_type` takes precedence.

- fast_sgd:

  If `TRUE`, then the following combination of parameters is set:
  `pcg_rand = TRUE`, `n_sgd_threads = "auto"` and `approx_pow = TRUE`.
  The default is `FALSE`. Setting this to `TRUE` will speed up the
  stochastic optimization phase, but give a potentially less accurate
  embedding, and which will not be exactly reproducible even with a
  fixed seed. For visualization, `fast_sgd = TRUE` will give perfectly
  good results. For more generic dimensionality reduction, it's safer to
  leave `fast_sgd = FALSE`. If `fast_sgd = TRUE`, then user-supplied
  values of `pcg_rand`, `n_sgd_threads`, and `approx_pow` are ignored.

- n_sgd_threads:

  Number of threads to use during stochastic gradient descent. If set to
  \> 1, then be aware that if `batch = FALSE`, results will *not* be
  reproducible, even if `set.seed` is called with a fixed seed before
  running. If set to `"auto"` then half the number of concurrent threads
  supported by the system will be used.

- grain_size:

  The minimum amount of work to do on each thread. If this value is set
  high enough, then less than `n_threads` or `n_sgd_threads` will be
  used for processing, which might give a performance improvement if the
  overhead of thread management and context switching was outweighing
  the improvement due to concurrent processing. This should be left at
  default (`1`) and work will be spread evenly over all the threads
  specified.

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
  (0/1) rather than real valued.

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

A matrix of optimized coordinates.

## References

Kingma, D. P., & Ba, J. (2014). Adam: A method for stochastic
optimization. *arXiv preprint* *arXiv*:1412.6980.
<https://arxiv.org/abs/1412.6980>

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

## Examples

``` r
iris30 <- iris[c(1:10, 51:60, 101:110), ]

# return a 30 x 30 sparse matrix with similarity data based on 10 nearest
# neighbors per item
iris30_sim_graph <- similarity_graph(iris30, n_neighbors = 10)
# produce 2D coordinates replicating the neighbor relations in the similarity
# graph
set.seed(42)
iris30_opt <- optimize_graph_layout(iris30_sim_graph, X = iris30)

# the above two steps are the same as:
# set.seed(42); iris_umap <- umap(iris30, n_neighbors = 10)
```
