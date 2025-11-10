# Add New Points to an Existing Embedding

Carry out an embedding of new data using an existing embedding. Requires
using the result of calling
[`umap`](https://jlmelville.github.io/uwot/reference/umap.md) or
[`tumap`](https://jlmelville.github.io/uwot/reference/tumap.md) with
`ret_model = TRUE`.

## Usage

``` r
umap_transform(
  X = NULL,
  model = NULL,
  nn_method = NULL,
  init_weighted = TRUE,
  search_k = NULL,
  tmpdir = tempdir(),
  n_epochs = NULL,
  n_threads = NULL,
  n_sgd_threads = 0,
  grain_size = 1,
  verbose = FALSE,
  init = "weighted",
  batch = NULL,
  learning_rate = NULL,
  opt_args = NULL,
  epoch_callback = NULL,
  ret_extra = NULL,
  seed = NULL
)
```

## Arguments

- X:

  The new data to be transformed, either a matrix of data frame. Must
  have the same columns in the same order as the input data used to
  generate the `model`.

- model:

  Data associated with an existing embedding.

- nn_method:

  Optional pre-calculated nearest neighbor data. There are two supported
  formats. The first is a list consisting of two elements:

  - `"idx"`. A `n_vertices x n_neighbors` matrix where `n_vertices` is
    the number of observations in `X`. The contents of the matrix should
    be the integer indexes of the data used to generate the `model`,
    which are the `n_neighbors`-nearest neighbors of the data to be
    transformed.

  - `"dist"`. A `n_vertices x n_neighbors` matrix containing the
    distances of the nearest neighbors.

  The second supported format is a sparse distance matrix of type
  `dgCMatrix`, with dimensions `n_model_vertices x n_vertices`. where
  `n_model_vertices` is the number of observations in the original data
  that generated the model. Distances should be arranged by column, i.e.
  a non-zero entry in row `j` of the `i`th column indicates that the
  `j`th observation in the original data used to generate the `model` is
  a nearest neighbor of the `i`th observation in the new data, with the
  distance given by the value of that element. In this format, a
  different number of neighbors is allowed for each observation, i.e.
  each column can contain a different number of non-zero values.
  Multiple nearest neighbor data (e.g. from two different pre-calculated
  metrics) can be passed by passing a list containing the nearest
  neighbor data lists as items.

- init_weighted:

  If `TRUE`, then initialize the embedded coordinates of `X` using a
  weighted average of the coordinates of the nearest neighbors from the
  original embedding in `model`, where the weights used are the edge
  weights from the UMAP smoothed knn distances. Otherwise, use an
  un-weighted average. This parameter will be deprecated and removed at
  version 1.0 of this package. Use the `init` parameter as a
  replacement, replacing `init_weighted = TRUE` with `init = "weighted"`
  and `init_weighted = FALSE` with `init = "average"`.

- search_k:

  Number of nodes to search during the neighbor retrieval. The larger k,
  the more the accurate results, but the longer the search takes.
  Default is the value used in building the `model` is used.

- tmpdir:

  Temporary directory to store nearest neighbor indexes during nearest
  neighbor search. Default is
  [`tempdir`](https://rdrr.io/r/base/tempfile.html). The index is only
  written to disk if `n_threads > 1`; otherwise, this parameter is
  ignored.

- n_epochs:

  Number of epochs to use during the optimization of the embedded
  coordinates. A value between `30 - 100` is a reasonable trade off
  between speed and thoroughness. By default, this value is set to one
  third the number of epochs used to build the `model`.

- n_threads:

  Number of threads to use, (except during stochastic gradient descent).
  Default is half the number of concurrent threads supported by the
  system.

- n_sgd_threads:

  Number of threads to use during stochastic gradient descent. If set to
  \> 1, then be aware that if `batch = FALSE`, results will *not* be
  reproducible, even if `set.seed` is called with a fixed seed before
  running. Set to `"auto"` to use the same value as `n_threads`.

- grain_size:

  Minimum batch size for multithreading. If the number of items to
  process in a thread falls below this number, then no threads will be
  used. Used in conjunction with `n_threads` and `n_sgd_threads`.

- verbose:

  If `TRUE`, log details to the console.

- init:

  how to initialize the transformed coordinates. One of:

  - `"weighted"` (The default). Use a weighted average of the
    coordinates of the nearest neighbors from the original embedding in
    `model`, where the weights used are the edge weights from the UMAP
    smoothed knn distances. Equivalent to `init_weighted = TRUE`.

  - `"average"`. Use the mean average of the coordinates of the nearest
    neighbors from the original embedding in `model`. Equivalent to
    `init_weighted = FALSE`.

  - A matrix of user-specified input coordinates, which must have
    dimensions the same as `(nrow(X), ncol(model$embedding))`.

  This parameter should be used in preference to `init_weighted`.

- batch:

  If `TRUE`, then embedding coordinates are updated at the end of each
  epoch rather than during the epoch. In batch mode, results are
  reproducible with a fixed random seed even with `n_sgd_threads > 1`,
  at the cost of a slightly higher memory use. You may also have to
  modify `learning_rate` and increase `n_epochs`, so whether this
  provides a speed increase over the single-threaded optimization is
  likely to be dataset and hardware-dependent. If `NULL`, the transform
  will use the value provided in the `model`, if available. Default:
  `FALSE`.

- learning_rate:

  Initial learning rate used in optimization of the coordinates. This
  overrides the value associated with the `model`. This should be left
  unspecified under most circumstances.

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

  If `NULL`, the transform will use the value provided in the `model`,
  if available.

- epoch_callback:

  A function which will be invoked at the end of every epoch. Its
  signature should be: `(epoch, n_epochs, coords, fixed_coords)`, where:

  - `epoch` The current epoch number (between `1` and `n_epochs`).

  - `n_epochs` Number of epochs to use during the optimization of the
    embedded coordinates.

  - `coords` The embedded coordinates as of the end of the current
    epoch, as a matrix with dimensions (N, `n_components`).

  - `fixed_coords` The originally embedded coordinates from the `model`.
    These are fixed and do not change. A matrix with dimensions (Nmodel,
    `n_components`) where `Nmodel` is the number of observations in the
    original data.

- ret_extra:

  A vector indicating what extra data to return. May contain any
  combination of the following strings:

  - `"fgraph"` the high dimensional fuzzy graph (i.e. the fuzzy
    simplicial set of the merged local views of the input data). The
    graph is returned as a sparse matrix of class
    [dgCMatrix-class](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)
    with dimensions `NX` x `Nmodel`, where `NX` is the number of items
    in the data to transform in `X`, and `NModel` is the number of items
    in the data used to build the UMAP `model`. A non-zero entry (i, j)
    gives the membership strength of the edge connecting the vertex
    representing the ith item in `X` to the jth item in the data used to
    build the `model`. Note that the graph is further sparsified by
    removing edges with sufficiently low membership strength that they
    would not be sampled by the probabilistic edge sampling employed for
    optimization and therefore the number of non-zero elements in the
    matrix is dependent on `n_epochs`. If you are only interested in the
    fuzzy input graph (e.g. for clustering), setting `n_epochs = 0` will
    avoid any further sparsifying.

  - `"nn"` the nearest neighbor graph for `X` with respect to the
    observations in the `model`. The graph will be returned as a list of
    two items: `idx` a matrix of indices, with as many rows as there are
    items in `X` and as many columns as there are nearest neighbors to
    be computed (this value is determined by the `model`). The indices
    are those of the rows of the data used to build the `model`, so
    they're not necessarily of much use unless you have access to that
    data. The second item, `dist` is a matrix of the equivalent
    distances, with the same dimensions as `idx`.

- seed:

  Integer seed to use to initialize the random number generator state.
  Combined with `n_sgd_threads = 1` or `batch = TRUE`, this should give
  consistent output across multiple runs on a given installation.
  Setting this value is equivalent to calling
  [`set.seed`](https://rdrr.io/r/base/Random.html), but it may be more
  convenient in some situations than having to call a separate function.
  The default is to not set a seed, in which case this function uses the
  behavior specified by the supplied `model`: If the model specifies a
  seed, then the model seed will be used to seed then random number
  generator, and results will still be consistent (if
  `n_sgd_threads = 1`). If you want to force the seed to not be set,
  even if it is set in `model`, set `seed = FALSE`.

## Value

A matrix of coordinates for `X` transformed into the space of the
`model`, or if `ret_extra` is specified, a list containing:

- `embedding` the matrix of optimized coordinates.

- if `ret_extra` contains `"fgraph"`, an item of the same name
  containing the high-dimensional fuzzy graph as a sparse matrix, of
  type
  [dgCMatrix-class](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html).

- if `ret_extra` contains `"sigma"`, returns a vector of the smooth knn
  distance normalization terms for each observation as `"sigma"` and a
  vector `"rho"` containing the largest distance to the locally
  connected neighbors of each observation.

- if `ret_extra` contains `"localr"`, an item of the same name
  containing a vector of the estimated local radii, the sum of `"sigma"`
  and `"rho"`.

- if `ret_extra` contains `"nn"`, an item of the same name containing
  the nearest neighbors of each item in `X` (with respect to the items
  that created the `model`).

## Details

Note that some settings are incompatible with the production of a UMAP
model via [`umap`](https://jlmelville.github.io/uwot/reference/umap.md):
external neighbor data (passed via a list to the argument of the
`nn_method` parameter), and factor columns that were included in the
UMAP calculation via the `metric` parameter. In the latter case, the
model produced is based only on the numeric data. A transformation is
possible, but factor columns in the new data are ignored.

## Examples

``` r
iris_train <- iris[1:100, ]
iris_test <- iris[101:150, ]

# You must set ret_model = TRUE to return extra data needed
iris_train_umap <- umap(iris_train, ret_model = TRUE)
iris_test_umap <- umap_transform(iris_test, iris_train_umap)
```
