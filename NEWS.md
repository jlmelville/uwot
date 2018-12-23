# uwot 0.0.0.9008 (December 23 2018)

## New features

* New parameter: `n_sgd_threads`, which controls the number of threads used
in the stochastic gradient descent. By default this is now single-threaded
and should result in reproducible results when using `set.seed`. To get back
the old, less consistent, but faster settings, set `n_sgd_threads = "auto"`.
* API change for consistency with Python UMAP:
  * `alpha` is now `learning_rate`.
  * `gamma` is now `repulsion_strength`.
* Default spectral initialization now looks for disconnected components and
initializes them separately (also applies to `laplacian` and `normlaplacian`).
* New `init` options: `sspectral`, `snormlaplacian` and `slaplacian`. These are
like `spectral`, `normlaplacian`, `laplacian` respectively, but scaled so that
each dimension has a standard deviation of 1e-4. This is like the difference
between the `pca` and `spca` options.

## Bug fixes and minor improvements

* Hamming distance support (was actually using Euclidean distance).
* Smooth knn/perplexity calibration results had a small dependency on the 
number of threads used.
* Anomalously long spectral intialization times should now be reduced.
* Internal changes and fixes thanks to a code review by Aaron Lun 
(https://github.com/ltla).

# uwot 0.0.0.9007 (December 9 2018)

## New features

* New parameter `pca`: set this to a positive integer to reduce matrix of
data frames to that number of columns using PCA. Only works if 
`metric = "euclidean"`. If you have > 100 columns, this can substantially 
improve the speed of the nearest neighbor search. t-SNE implementations often
set this value to 50.

## Bug fixes and minor improvements

* Laplacian Eigenmap initialization convergence failure is now correctly 
detected.
* C++ code was over-writing data passed from R as a function argument.

# uwot 0.0.0.9006 (December 5 2018)

## New features

* Highly experimental mixed data type support for `metric`: instead of
specifying a single metric name (e.g. `metric = "euclidean"`), you can pass a
list, where the name of each item is the metric to use and the value is a vector
of the names of the columns to use with that metric, e.g. 
`metric = list("euclidean" = c("A1", "A2"), "cosine" = c("B1", "B2", "B3"))` 
treats columns `A1` and `A2` as one block, using the Euclidean distance to find
nearest neighbors, whereas `B1`, `B2` and `B3` are treated as a second block,
using the cosine distance.
* Factor columns can also be used in the metric, using the metric name 
`categorical`. 
* `y` may now be a data frame or matrix if multiple target data is available. 
* New parameter `target_metric`, to specify the distance metric to use with 
numerical `y`. This has the same capabilities as `metric`.
* Multiple external nearest neighbor data sources are now supported. Instead of
passing a list of two matrices, pass a list of lists, one for each external
metric.
* More details on mixed data types can be found at 
https://github.com/jlmelville/uwot#mixed-data-types.
* Compatibility with older versions of RcppParallel (contributed by 
[sirusb](https://github.com/sirusb)).
* `scale = "Z"` To Z-scale each column of input (synonym for `scale = TRUE` 
or `scale = "scale"`).
* New scaling option, `scale = "colrange"` to scale columns in the range (0, 1).

# uwot 0.0.0.9005 (November 4 2018)

## New features

* Hamming distance is now supported, due to upgrade to RcppAnnoy 0.0.11.

# uwot 0.0.0.9004 (October 21 2018)

## New features

* For supervised UMAP with numeric `y`, you may pass nearest neighbor data
directly, in the same format as that supported by `X`-related nearest neighbor
data. This may be useful if you don't want to use Euclidean distances for
the `y` data, or if you have missing data (and have a way to assign nearest neighbors
for those cases, obviously). See the 
[Nearest Neighbor Data Format](https://github.com/jlmelville/uwot#nearest-neighbor-data-format)
section for details.

# uwot 0.0.0.9003 (September 22 2018)

## New features

* New parameter `ret_nn`: when `TRUE` returns nearest neighbor matrices
as a `nn` list: indices in item `idx` and distances in item `dist`. Embedded
coordinates are in `embedding`. Both `ret_nn` and `ret_model` can be `TRUE`,
and should not cause any compatibility issues with supervised embeddings.
* `nn_method` can now take precomputed nearest neighbor data. Must be a list of
two matrices: `idx`, containing integer indexes, and `dist` containing 
distances. By no coincidence, this is the format return by `ret_nn`.

## Bug fixes and minor improvements

* Embedding to `n_components = 1` was broken (https://github.com/jlmelville/uwot/issues/6)
* User-supplied matrices to `init` parameter were being modified, in defiance of basic R pass-by-copy semantics.

# uwot 0.0.0.9002 (August 14 2018)

## Bug fixes and minor improvements

* `metric = "cosine"` is working again for `n_threads` greater than `0` (https://github.com/jlmelville/uwot/issues/5)

# uwot 0.0.0.9001

## New features

* *August 5 2018*. You can now use an existing embedding to add new points via
`umap_transform`. See the example section below.

* *August 1 2018*. Numerical vectors are now supported for supervised dimension reduction.

* *July 31 2018*. (Very) initial support for supervised dimension reduction:
categorical data only at the moment. Pass in a factor vector (use `NA` for
unknown labels) as the `y` parameter and edges with bad (or unknown) labels are
down-weighted, hopefully leading to better separation of classes. This works
remarkably well for the Fashion MNIST dataset.

* *July 22 2018*. You can now use the cosine and Manhattan distances with the
Annoy nearest neighbor search, via `metric = "cosine"` and `metric =
"manhattan"`, respectively. Hamming distance is not supported because RcppAnnoy
doesn't yet support it.
