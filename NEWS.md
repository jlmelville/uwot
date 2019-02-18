# uwot 0.0.0.9010

## New features

* New parameter: `n_refine_iters`, which specifies how many iterations of
[nearest neighbor descent](https://doi.org/10.1145/1963405.1963487) to carry out
to refine the results of the Annoy approximate neighbor search. When set to 
greater than zero, the `search_k` parameter is reduced. Default is currently 0,
so it's not currently on by default, because although the accuracy of the
neighbor data is likely to be improved, the overall run time is slightly
increased. If you want to use it, LargeVis uses 3 iterations, which is a good
place to start.
* New parameter: `init_sdev` which specifies how large the standard deviation
of each column of the initial coordinates should be. This will scale any input
coordinates (including user-provided matrix coordinates). `init = "spca"` can
now be thought of as an alias of `init = "pca", init_sdev = 1e-4`. 
This may be too aggressive scaling for some datasets. The typical UMAP spectral
initializations tend to result in standard deviations of around `2` to `5`, so
this might be more appropriate in some cases.
* As a result of adding `init_sdev`, the `init` options `sspectral`,
`slaplacian` and `snormlaplacian` have been removed (they weren't around for
very long anyway). You can get the same behavior by e.g. 
`init = "spectral", init_sdev = 1e-4`. `init = "spca"` is sticking around 
because I use it a lot.

## Bug fixes and minor improvements

* If requesting a spectral initialization, but multiple disconnected components
are present, fall back to `init = "spca"`. 
* Removed dependency on C++ `<random>` header. This breaks backwards
compatibility in the sense that results with a given seed will not be the same
as with previous versions.
* `metric = "cosine"` results were incorrectly using the unmodified Annoy
angular distance.
* Numeric matrix columns can be specified as the target for the `categorical`
metric (fixes https://github.com/jlmelville/uwot/issues/20).

# uwot 0.0.0.9009 (1 January 2019)

* Data is now stored column-wise during optimization, which should result in
an increase in performance for larger values of `n_components` (e.g. 
approximately 50% faster optimization time with MNIST and `n_components = 50`).
* New parameter: `pca_center`, which controls whether to center the data before
applying PCA. It would be typical to set this to `FALSE` if you are applying
PCA to binary data (although note you can't use this with setting with 
`metric = "hamming"`)
* PCA will now be used when the `metric` is `"manhattan"` and `"cosine"`. It's
still *not* applied when using `"hamming"` (data still needs to be in binary
format, not real-valued).
* If using mixed datatypes, you may override the `pca` and `pca_center`
parameter values for a given data block by using a list for the value of the
metric, with the column ids/names as an unnamed item and the overriding values
as named items, e.g. instead of `manhattan = 1:100`, use 
`manhattan = list(1:100, pca_center = FALSE)` to turn off PCA centering for 
just that block. This functionality exists mainly for the case where you have
mixed binary and real-valued data and want to apply PCA to both data types. It's
normal to apply centering to real-valued data but not to binary data.

## Bug fixes and minor improvements

* Fixed bug that affected `umap_transform`, where negative sampling was over 
the size of the test data (should be the training data).
* Some other performance improvements (around 10% faster for the optimization
stage with MNIST).
* When `verbose = TRUE`, log the Annoy recall accuracy, which may help tune
values of `n_trees` and `search_k`.

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
* Anomalously long spectral initialization times should now be reduced.
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
