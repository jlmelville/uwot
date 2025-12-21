# uwot (development version)

## Bug fixes and minor improvements

* New parameter: `n_build_threads`. Controls the number of threads used to build
the nearest neighbor index when using `nn_method = "hnsw"` or
`nn_method = "nndescent"`. With a given seed, setting `n_build_threads = 1` will
give deterministic nearest neighbor results when using either of those methods
(at the cost of some speed). If you are using the default Annoy-based method,
this parameter does nothing: it has always used one thread to build the
index. Thank you [hsuknowledge](https://github.com/hsuknowledge) for pointing
out the issue around determinism.
(<https://github.com/jlmelville/uwot/issues/139>).
* When passing a list of two matrices to `nn_method` to represent a pre-computed
k-nearest neighbor graph, you may specify the indices and distances as `index`
and `distance` respectively. This improves compatibility with the
`BiocNeighbors` package. Thank you [hsuknowledge](https://github.com/hsuknowledge)
for pointing out the mismatch in naming convention
(<https://github.com/jlmelville/uwot/issues/140>).

# uwot 0.2.4

## Bug fixes and minor improvements
 
* The installation status of optional dependencies were not being detected
correctly. This meant that different packages could be used for initialization
in unpredictable ways depending on whether they had been explicitly loaded or
not. Thank you [hsuknowledge](https://github.com/hsuknowledge) for the report
(<https://github.com/jlmelville/uwot/issues/134>).
* Users of the [bbknnR](https://github.com/ycli1995/bbknnR/issues/3) package 
were suffering from a not-helpful error message when the custom neighbor data
contained missing values. An explicit check has been added and although this is
still a fatal error, the message should be more informative
(<https://github.com/jlmelville/uwot/issues/135>).
* Fixed a partially-specified parameter name being passed to irlba. Thank you
[Hugo Gruson](https://github.com/Bisaloo) for the fix 
(<https://github.com/jlmelville/uwot/pull/136>).
* Development-only: fixed an incorrect use of `testthat::expect` in a unit test.
Thank you [Hadley Wickham](https://github.com/hadley) for the fix
(<https://github.com/jlmelville/uwot/pull/138>).

# uwot 0.2.3

## New features:

* New parameter: `rng_type`. This will be used in favor of the boolean 
`pcg_rand` parameter, although `pcg_rand` will still work for backwards
compatibility.
* New negative sampling option: set `rng_type = "deterministic"` to use a
deterministic sampling of vertices during the optimization phase. This should
give qualitatively similar results to using a real PRNG, but has the advantage
of being faster and giving more reproducible output. This feature was inspired
by a comment by 
[Leland McInnes on Reddit](https://www.reddit.com/r/MachineLearning/comments/1gsjfq9/comment/lxip9wy/).

## Bug fixes and minor improvements

* Setting `num_threads` directly in `umap2` did not result in the number of SGD
threads being updated to that value when `batch = TRUE`, which it should have
been.
* Despite assertions to the contrary in version 0.2.1, `umap_transform`
continued to return the fuzzy graph in transposed form. Thank you
[PedroMilanezAlmeida](https://github.com/PedroMilanezAlmeida) for 
reopening the issue (<https://github.com/jlmelville/uwot/issues/118>).
* Relative paths could not be used to save a model. Thank you
[Wouter van der Bijl](https://github.com/Ax3man) for the bug report
(<https://github.com/jlmelville/uwot/issues/131>) and the suggested fix.
* `repulsion_strength` was silently ignored if used with `tumap` or `umap2` with
`a = 1, b = 1`. Ignoring the setting was on purpose, but it was not documented
anywhere. `repulsion_strength` is now compatible with these settings.
* It's no longer an error to provide a `pca` argument if the input data has a
maximum rank smaller than the value of `pca`. No PCA is applied in this case.
If `verbose = TRUE`, a message will be printed to inform the user.

# uwot 0.2.2

## Bug fixes and minor improvements

* `RSpectra` is now a required dependency (again). It was a required dependency
up until version 0.1.12, when it became optional (`irlba` was used in its
place). However, problems with interactions of the current version of `irlba`
with an ABI change in the `Matrix` package means that it's hard for downstream
packages and users to build `uwot` without re-installing `Matrix` and `irlba`
from source, which may not be an option for some people. Also it was causing a
CRAN check error. I have changed some tests, examples and vignettes to use
`RSpectra` explicitly, and to only test `irlba` code-paths where necessary. See
<https://github.com/jlmelville/uwot/issues/115> and links therein for more
details.

# uwot 0.2.1

## New features:

* The [HNSW](https://github.com/nmslib/hnswlib) approximate nearest neighbor
search algorithm is now supported via the
[RcppHNSW](https://cran.r-project.org/package=RcppHNSW) package. Set
`nn_method = "hnsw"` to use it. The behavior of the method can be controlled by
the new `nn_args` parameter, a list which may contain `M`, `ef_construction`
and `ef`. See the hnswlib library's 
[ALGO_PARAMS documentation](https://github.com/nmslib/hnswlib/blob/master/ALGO_PARAMS.md)
for details on these parameters. Although typically faster than Annoy (for a
given accuracy), be aware that the only supported `metric` values are
`"euclidean"`, `"cosine"` and `"correlation"`. Finally, RcppHNSW is only a
suggested package, not a requirement, so you need to install it yourself (e.g.
via `install.packages("RcppHNSW")`). Also see the
[article on HNSW in uwot](https://jlmelville.github.io/uwot/articles/hnsw-umap.html)
in the documentation.
* The nearest neighbor descent approximate nearest neighbor search algorithm is
now supported via the
[rnndescent](https://cran.r-project.org/package=rnndescent) package. Set
`nn_method = "nndescent"` to use it. The behavior of the method can be
controlled by the new `nn_args` parameter. There are many supported metrics and
possible parameters that can be set in `nn_args`, so please see the 
[article on nearest neighbor descent in uwot](https://jlmelville.github.io/uwot/articles/rnndescent-umap.html)
in the documentation, and also the rnndescent package's
[documentation](https://jlmelville.github.io/rnndescent/index.html) for details.
`rnndescent` is only a suggested package, not a requirement, so you need to 
install it yourself (e.g. via `install.packages("rnndescent")`).
* New function: `umap2`, which acts like `umap` but with modified defaults,
reflecting my experience with UMAP and correcting some small mistakes. See the
[umap2 article](https://jlmelville.github.io/uwot/articles/umap2.html) for more
details.

## Bug fixes and minor improvements

* `init_sdev = "range"` caused an error with a user-supplied `init` matrix.
* Transforming new data with the `correlation` metric was actually using the
`cosine` metric if you saved and reloaded the model. Thank you 
[Holly Hall](https://github.com/mdrnao) for the report and helpful detective
work (<https://github.com/jlmelville/uwot/issues/117>).
* `umap_transform` could fail if the new data to be transformed had the 
`scaled:center` and `scaled:scale` attributes set (e.g. from applying the
`scale` function).
* If you asked `umap_transform` to return the fuzzy graph (
`ret_extra = c("fgraph")`), it was transposed when `batch = TRUE, n_epochs = 0`.
Thank you [PedroMilanezAlmeida](https://github.com/PedroMilanezAlmeida) for 
reporting (<https://github.com/jlmelville/uwot/issues/118>).
* Setting `n_sgd_threads = "auto"` with `umap_transform` caused a crash.
* A warning was being emitted due to not being specific enough about what `dist`
class was meant that may have been particularly affecting Seurat users. Thank 
you [AndiMunteanu](https://github.com/AndiMunteanu) for reporting (and
suggesting a solution) (<https://github.com/jlmelville/uwot/issues/121>).

# uwot 0.1.16

## Bug fixes and minor improvements

* A small change to a header file was required to fully support the next version
of [RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy). Thank
you [Dirk Eddelbuettel](https://github.com/eddelbuettel) for the PR
(<https://github.com/jlmelville/uwot/issues/112>).

# uwot 0.1.15

## New features:

* New function: `optimize_graph_layout`. Use this to produce optimized output 
coordinates that reflect an input similarity graph (such as that produced by
the `similarity_graph` function. `similarity_graph` followed by 
`optimize_graph_layout` is the same as running `umap`, so the purpose of these
functions is to allow for more flexibility and decoupling between generating
the nearest neighbor graph and optimizing the low-dimensional approximation
to it. Based on a request by user [Chengwei94](https://github.com/Chengwei94) 
(<https://github.com/jlmelville/uwot/issues/98>).
* New functions: `simplicial_set_union` and `simplicial_set_intersect`. These
allow for the combination of different fuzzy graph representations of a dataset
into a single fuzzy graph using the UMAP simplicial set operations. Based on a
request in the Python UMAP issues tracker by user 
[Dhar xion](https://github.com/ratheraarif).
* New parameter for `umap_transform`: `ret_extra`. This works like the
equivalent parameter for `umap`, and should be a character vector specifying the
extra information you would like returned in addition to the embedding, in which
case a list will be returned with an `embedding` member containing the optimized
coordinates. Supported values are `"fgraph"`, `"nn"`, `"sigma"` and `"localr"`.
Based on a request by user 
[PedroMilanezAlmeida](https://github.com/PedroMilanezAlmeida) 
(<https://github.com/jlmelville/uwot/issues/104>).
* New parameter from `umap`, `tumap` and `umap_transform`: `seed`. This will do
the equivalent of calling `set.seed` internally, and hence will help with 
reproducibility. The chosen seed is exported if `ret_model = TRUE` and 
`umap_transform` will use that seed if present, so you only need to specify
it in `umap_transform` if you want to change the seed. The default behavior 
remains to not modify the random number state. Based on a request by
[SuhasSrinivasan](https://github.com/SuhasSrinivasan)
(<https://github.com/jlmelville/uwot/issues/110>).

## Bug fixes and minor improvements

* A new setting for `init_sdev`: set `init_sdev = "range"` and initial
coordinates will be range-scaled so each column takes values between 0-10. This
pre-processing was added to the Python UMAP package at some point after `uwot`
began development and so should probably always be used with the default 
`init = "spectral"` setting. However, it is not set by default to maintain
backwards compatibility with older versions of `uwot`.
* `ret_extra = c("sigma")` is now supported by `lvish`. The Gaussian bandwidths
are returned in a `sigma` vector. In addition, a vector of intrinsic
dimensionalities estimated for each point using an analytical expression of the
finite difference method given by 
[Lee and co-workers](https://doi.org/10.1016/j.neucom.2014.12.095) is returned 
in the `dint` vector.
* The `min_dist` and `spread` parameters are now returned in the model when
`umap` is run with `ret_model = TRUE`. This is just for documentation purposes, 
these values are not used directly by the model in `umap_transform`. If the 
parameters `a` and `b` are set directly when invoking `umap`, then both 
`min_dist` and `spread` will be set to `NULL` in the returned model. This 
feature was added in response to a question from 
[kjiang18](https://github.com/kjiang18) 
(<https://github.com/jlmelville/uwot/issues/95>).
* Some new checks for NA values in input data have been added. Also a warning
will be emitted if `n_components` seems to have been set too high.
* If `n_components` was greater than `n_neighbors` then `umap_transform` would
crash the R session. Thank you to [ChVav](https://github.com/ChVav) for 
reporting this (<https://github.com/jlmelville/uwot/issues/102>).
* Using `umap_transform` with a model where `dens_scale` was set could cause
a segmentation fault, destroying the session. Even if it didn't it could give
an entirely artifactual "ring" structure. Thank you 
[FemkeSmit](https://github.com/FemkeSmit) for reporting this and providing
assistance in diagnosing the underlying cause 
(<https://github.com/jlmelville/uwot/issues/103>).
* If you set `binary_edge_weights = TRUE`, this setting was not exported when
`ret_model = TRUE`, and was therefore not respected by `umap_transform`. This 
has now been fixed, but you will need to regenerate any models that used
binary edge weights.
* The rdoc for the `init` param said that if there were multiple disconnected
components, a spectral initialization would attempt to merge multiple 
sub-graphs. Not true: actually, spectral initialization is abandoned in favor
of PCA. The documentation has been updated to reflect the true state of affairs.
No idea what I was thinking of there.
* `load_model` and `save_model` didn't work on Windows 7 due to how the version
of `tar` there handles drive letters. Thank you
[mytarmail](https://github.com/mytarmail) for the report 
(<https://github.com/jlmelville/uwot/issues/109>).
* Warn if the initial coordinates have a very large scale (a standard deviation 
> 10.0), because this can lead to small gradients and poor optimization. Thank
you [SuhasSrinivasan](https://github.com/SuhasSrinivasan) for the report
(<https://github.com/jlmelville/uwot/issues/110>).
* A change to accommodate a forthcoming version of 
[RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy). Thank
you [Dirk Eddelbuettel](https://github.com/eddelbuettel) for the PR
(<https://github.com/jlmelville/uwot/issues/111>).

# uwot 0.1.14

## New features

* New function: `similarity_graph`. If you are more interested in the
high-dimensional graph/fuzzy simplicial set representation of your input data,
and don't care about the low dimensional approximation, the `similarity_graph`
function offers a similar API to `umap`, but neither the initialization nor
optimization of low-dimensional coordinates will be performed. The return value
is the same as that which would be returned in the results list as the `fgraph`
member if you had provided `ret_extra = c("fgraph")`. Compared to getting the
same result via running `umap`, this function is a bit more convenient to use,
makes your intention clearer if you would be discarding the embedding, and saves
a small amount of time. A t-SNE/LargeVis similarity graph can be returned by
setting `method = "largevis"`.

## Bug fixes and minor improvements

* If a model was generated without using pre-generated nearest neighbors, you
couldn't use `umap_transform` with pre-generated nearest neighbors (also the
error message was completely useless). Thank you to 
[AustinHartman](https://github.com/AustinHartman) for reporting this 
(<https://github.com/jlmelville/uwot/issues/97>).

# uwot 0.1.13

* This is a resubmission of 0.1.12 but with an internal function 
(`fuzzy_simplicial_set`) refactored to behave more like that of previous 
versions. This change was breaking the behavior of the CRAN package 
[bbknnR](https://cran.r-project.org/package=bbknnR).

# uwot 0.1.12

## New features

* New parameter: `dens_weight`. If set to a value between 0 and 1, an attempt
is made to include the relative local densities of the input data in the output
coordinates. This is an approximation to the 
[densMAP](https://doi.org/10.1038/s41587-020-00801-7) method. A large value of
`dens_weight` will use a larger range of output densities to reflect the input
data. If the data is too spread out, reduce the value of `dens_weight`. For
more information see the 
[documentation at the uwot repo](https://jlmelville.github.io/uwot/articles/leopold.html).
* New parameter: `binary_edge_weights`. If set to `TRUE`, instead of smoothed
knn distances, non-zero edge weights all have a value of 1. This is how
[PaCMAP](https://www.jmlr.org/papers/v22/20-1061.html) works and there is
[practical](https://arxiv.org/abs/2007.08902) and
[theoretical](https://proceedings.neurips.cc/paper/2021/hash/2de5d16682c3c35007e4e92982f1a2ba-Abstract.html)
reasons to believe this won't have a big effect on UMAP but you can try it
yourself.
* New options for `ret_extra`: 
    * `"sigma"`: the return value will contain a `sigma` entry, a vector of the
    smooth knn distance scaling normalization factors, one for each observation
    in the input data. A small value indicates a high density of points in the
    local neighborhood of that observation. For `lvish` the equivalent
    bandwidths calculated for the input perplexity is returned.
    * also, a vector `rho` will be exported, which is the distance to the
    nearest neighbor after the number of neighbors specified by the
    `local_connectivity`. Only applies for `umap` and `tumap`.
    * `"localr"`: exports a vector of the local radii, the sum of `sigma` and
    `rho` and used to scale the output coordinates when `dens_weight` is set.
    Even if not using `dens_weight`, visualizing the output coordinates using a
    color scale based on the value of `localr` can reveal regions of the input
    data with different densities.
* For functions `umap` and `tumap` only: new data type for precomputed nearest
neighbor data passed as the `nn_method` parameter: you may use a sparse distance
matrix of format `dgCMatrix` with dimensions `N x N` where `N` is the number of
observations in the input data. Distances should be arranged by column, i.e. a
non-zero entry in row `j` of the `i`th column indicates that the `j`th
observation in the input data is a nearest neighbor of the `i`th observation
with the distance given by the value of that element. Note that this is a
different format to the sparse distance matrix that can be passed as input to
`X`: notably, the matrix is not assumed to be symmetric. Unlike other input
formats, you may have a different number of neighbors for each observation (but
there must be at least one neighbor defined per observation).
* `umap_transform` can also take a sparse distance matrix as its `nn_method`
parameter if precomputed nearest neighbor data is used to generate an initial
model. The format is the same as for the `nn_method` with `umap`. Because
distances are arranged by columns, the expected dimensions of the sparse matrix
is `N_model x N_new` where `N_model` is the number of observations in the
original data and `N_new` is the number of observations in the data to be
transformed.

## Bug fixes and minor improvements

* Models couldn't be re-saved after loading. Thank you to 
[ilyakorsunsky](https://github.com/ilyakorsunsky) for reporting this 
(<https://github.com/jlmelville/uwot/issues/88>).
* [RSpectra](https://cran.r-project.org/package=RSpectra) is now a 'Suggests',
rather than an 'Imports'. If you have RSpectra installed, it will be used
automatically where previous versions required it (for spectral initialization).
Otherwise, [irlba](https://cran.r-project.org/package=irlba) will be used.
For two-dimensional output, you are unlikely to notice much difference in speed
or accuracy with real-world data. For highly-structured simulation datasets
(e.g. spectral initialization of a 1D line) then RSpectra will give much better,
faster initializations, but these are not the typical use cases envisaged for
this package. For embedding into higher dimensions (e.g. `n_components = 100` or
higher), RSpectra is recommended and will likely out-perform irlba even if you
have installed a good linear algebra library.
* `init = "laplacian"` returned the wrong coordinates because of a slightly 
subtle issue around how to order the eigenvectors when using the random walk
transition matrix rather than normalized graph laplacians.
* The `init_sdev` parameter was ignored when the `init` parameter was a
user-supplied matrix. Now the input will be scaled.
* Matrix input was being converted to and from a data frame during
pre-processing, causing R to allocate memory that it was disinclined to ever
give up even after the function exited. This unnecessary manipulation is now
avoided.
* The behavior of the `bandwidth` parameter has been changed to give results
more like the current version (0.5.2) of the Python UMAP implementation. This is
likely to be a breaking change for non-default settings of `bandwidth`, but this
is not a parameter which is actually exposed by the Python UMAP public API any
more, so is on the road to deprecation in uwot too and I don't recommend you
change this.
* Transforming data with multiple blocks would give an error if the number of
rows of the new data did not equal the number of number of rows in the original
data.

# uwot 0.1.11

## New features

* New parameter: `batch`. If `TRUE`, then results are reproducible when
`n_sgd_threads > 1` (as long as you use `set.seed`). The price to be paid is
that the optimization is slightly less efficient (because coordinates are not
updated as quickly and hence gradients are staler for longer), so it is highly
recommended to set `n_epochs = 500` or higher. Thank you to 
[Aaron Lun](https://github.com/LTLA) who not only came up with a way to
implement this feature, but also wrote an entire 
[C++ implementation of UMAP](https://github.com/libscran/umappp) which does it 
(<https://github.com/jlmelville/uwot/issues/83>).
* New parameter: `opt_args`. The default optimization method when `batch = TRUE`
is [Adam](https://arxiv.org/abs/1412.6980). You can control its parameters by
passing them in the `opt_args` list. As Adam is a momentum-based method it
requires extra storage of previous gradient data. To avoid the extra memory
overhead you can also use `opt_args = list(method = "sgd")` to use a stochastic
gradient descent method like that used when `batch = FALSE`.
* New parameter: `epoch_callback`. You may now pass a function which will be
invoked at the end of each epoch. Mainly useful for producing an image of the
state of the embedding at different points during the optimization. This is
another feature taken from [umappp](https://github.com/libscran/umappp).
* New parameter: `pca_method`, used when the `pca` parameter is supplied to
reduce the initial dimensionality of the data. This controls which method is
used to carry out the PCA and can be set to one of:
    * `"irlba"` which uses `irlba::irlba` to calculate a truncated SVD. If this
    routine deems that you are trying to extract 50% or more of the singular 
    vectors, you will see a warning to that effect logged to the console.
    * `"rsvd"`, which uses `irlba::svdr` for truncated SVD. This method uses a
    small number of iterations which should give an accuracy/speed up trade-off
    similar to that of the 
    [scikit-learn TruncatedSVD](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.TruncatedSVD.html#sklearn.decomposition.TruncatedSVD)
    method. This can be much faster than using `"irlba"` but potentially at a
    cost in accuracy. However, for the purposes of dimensionality reduction as
    input to nearest neighbor search, this doesn't seem to matter much.
    * `"bigstatsr"`, which uses the [bigstatsr](https://cran.r-project.org/package=bigstatsr) 
    package will be used. **Note**: that this is *not* a dependency of `uwot`.
    If you want to use `bigstatsr`, you must install it yourself. On platforms
    without easy access to fast linear algebra libraries (e.g. Windows), using
    `bigstatsr` may give a speed up to PCA calculations.
    * `"svd"`, which uses `base::svd`. **Warning**: this is likely to be very
    slow for most datasets and exists as a fallback for small datasets where
    the `"irlba"` method would print a warning.
    * `"auto"` (the default) which uses `"irlba"` to calculate a truncated SVD,
    unless you are attempting to extract 50% or more of the singular vectors,
    in which case `"svd"` is used.

## Bug fixes and minor improvements

* If row names are provided in the input data (or nearest neighbor data, or
initialization data if it's a matrix), this will be used to name the rows of the
output embedding (<https://github.com/jlmelville/uwot/issues/81>), and also the
nearest neighbor data if you set `ret_nn = TRUE`. If the names exist in more
than one of the input data parameters listed above, but are inconsistent, no
guarantees are made about which names will be used. Thank you 
[jwijffels](https://github.com/jwijffels) for reporting this.
* In `umap_transform`, the learning rate is now down-scaled by a factor of 4,
consistent with the Python implementation of UMAP. If you need the old behavior
back, use the (newly added) `learning_rate` parameter in `umap_transform` to set
it explicitly. If you used the default value in `umap` when creating the model,
the correct setting in `umap_transform` is `learning_rate = 1.0`.
* Setting `nn_method = "annoy"` and `verbose = TRUE` would lead to an error with 
datasets with fewer than 50 items in them.
* Using multiple pre-computed nearest neighbors blocks is now supported with
`umap_transform` (this was incorrectly documented to work).
* Documentation around pre-calculated nearest neighbor data for `umap_transform`
was wrong in other ways: it has now been corrected to indicate that there should
be neighbor data for each item in the test data, but the neighbors and distances
should refer to items in training data (i.e. the data used to build the model).
* `n_neighbors` parameter is now correctly ignored in model generation if
pre-calculated nearest neighbor data is provided.
* Documentation incorrectly said `grain_size` didn't do anything.

# uwot 0.1.10

This release is mainly to allow for some internal changes to keep compatibility
with RcppAnnoy, used for the nearest neighbor calculations.

## Bug fixes and minor improvements

* Passing in data with missing values will now raise an error early. Missing
data in factor columns intended for supervised UMAP is still ok. Thank you David
McGaughey for tweeting about this issue.
* The documentation for the return value of `umap` and `tumap` now note that the
contents of the `model` list are subject to change and not intended to be part
of the uwot public API. I recommend not relying on the structure of the `model`,
especially if your package is intended to appear on CRAN or Bioconductor, as any
breakages will delay future releases of uwot to CRAN.

# uwot 0.1.9

## New features

* New metric: `metric = "correlation"` a distance based on the Pearson 
correlation (<https://github.com/jlmelville/uwot/issues/22>). Supporting this
required a change to the internals of how nearest neighbor data is stored.
Backwards compatibility with models generated by previous versions using 
`ret_model = TRUE` should have been preserved.

## Bug fixes and minor improvements

* New parameter, `nn_method`, for `umap_transform`: pass a list containing
pre-computed nearest neighbor data (identical to that used in the `umap`
function). You should not pass anything to the `X` parameter in this case. This
extends the functionality for transforming new points to the case where nearest
neighbor data between the original data and new data can be calculated external
to `uwot`. Thanks to [Yuhan Hao](https://github.com/yuhanH) for contributing the
PR (<https://github.com/jlmelville/uwot/issues/63> and
<https://github.com/jlmelville/uwot/issues/64>).
* New parameter, `init`, for `umap_transform`: provides a variety of options for
initializing the output coordinates, analogously to the same parameter in the
`umap` function (but without as many options currently). This is intended to
replace `init_weighted`, which should be considered deprecated, but won't be
removed until uwot 1.0 (whenever that is). Instead of `init_weighted = TRUE`,
use `init = "weighted"`; replace `init_weighted = FALSE` with 
`init = "average"`. Additionally, you can pass a matrix to `init` to act as the 
initial coordinates.
* Also in `umap_transform`: previously, setting `n_epochs = 0` was ignored: at
least one iteration of optimization was applied. Now, `n_epochs = 0` is
respected, and will return the initialized coordinates without any further
optimization.
* Minor performance improvement for single-threaded nearest neighbor search when
`verbose = TRUE`: the progress bar calculations were taking up a detectable
amount of time and has now been fixed. With very small data sets (< 50 items) the
progress bar will no longer appear when building the index.
* Passing a sparse distance matrix as input now supports upper/lower triangular
matrix storage rather than wasting storage using an explicitly symmetric 
sparse matrix.
* Minor license change: uwot used to be licensed under GPL-3 only; now it is 
GPL-3 or later.

# uwot 0.1.8

## Bug fixes and minor improvements

* default for `n_threads` is now `NULL` to provide a bit more protection from
changing dependencies.
* parallel code now uses the standard C++11 implementation of threading rather 
than tinythread++.
* The `grain_size` parameter has been undeprecated. As the version that 
deprecated this never made it to CRAN, this is unlikely to have affected many
people.

# uwot 0.1.7

## Bug fixes and minor improvements

* uwot should no longer trigger undefined behavior in sanitizers, due to the
temporary replacement of the RcppParallel package with code "borrowed" from
that package and using tinythread++ rather than tbb 
(<https://github.com/jlmelville/uwot/issues/52>).
* Further sanitizer improvements in the nearest neighbor search code due to
the upstream efforts of [erikbern](https://github.com/erikbern) and
[eddelbuettel](https://github.com/eddelbuettel)
(<https://github.com/jlmelville/uwot/issues/50>).
* The `grain_size` parameter is now ignored and remains to avoid breaking
backwards compatibility only.

# uwot 0.1.6

## New features

* New parameter, `ret_extra`, a vector which can contain any combination of: 
`"model"` (same as `ret_model = TRUE`), `"nn"` (same as `ret_nn = TRUE`) and 
`fgraph` (see below).
* New return value data: If the `ret_extra` vector contains `"fgraph"`, the 
returned list will contain an `fgraph` item representing the fuzzy simplicial 
input graph as a sparse N x N matrix. For `lvish`, use `"P"` instead of 
`"fgraph`" (<https://github.com/jlmelville/uwot/issues/47>). Note that there
is a further sparsifying step where edges with a very low membership are removed
if there is no prospect of the edge being sampled during optimization. This is
controlled by `n_epochs`: the smaller the value, the more sparsifying will
occur. If you are only interested in the fuzzy graph and not the embedded
coordinates, set `n_epochs = 0`.
* New function: `unload_uwot`, to unload the Annoy nearest neighbor indices in
a model. This prevents the model from being used in `umap_transform`, but allows
for the temporary working directory created by both `save_uwot` and `load_uwot`
to be deleted. Previously, both `load_uwot` and `save_uwot` were attempting to
delete the temporary working directories they used, but would always silently
fail because Annoy is making use of files in those directories.
* An attempt has been made to reduce the variability of results due to different
compiler and C++ library versions on different machines. Visually results are
unchanged in most cases, but this is a breaking change in terms of numerical
output. The best chance of obtaining floating point determinism across machines
is to use `init = "spca"`, fixed values of `a` and `b` (rather than allowing
them to be calculated through setting `min_dist` and `spread`) and
`approx_pow = TRUE`. Using the `tumap` method with `init = "spca"` is probably
the most robust approach.

## Bug fixes and minor improvements

* New behavior when `n_epochs = 0`. This used to behave like (`n_epochs = NULL`)
and gave a default number of epochs (dependent on the number of vertices in the
dataset). Now it more usefully carries out all calculations except optimization,
so the returned coordinates are those specified by the `init` parameter, so this
is an easy way to access e.g. the spectral or PCA initialization coordinates.
If you want the input fuzzy graph (`ret_extra` vector contains `"fgraph"`), this
will also prevent the graph having edges with very low membership being removed.
You still get the old default epochs behavior by setting `n_epochs = NULL` or to
a negative value.
* `save_uwot` and `load_uwot` have been updated with a `verbose` parameter so
it's easier to see what temporary files are being created.
* `save_uwot` has a new parameter, `unload`, which if set to `TRUE` will delete
the working directory for you, at the cost of unloading the model, i.e. it can't 
be used with `umap_transform` until you reload it with `load_uwot`.
* `save_uwot` now returns the saved model with an extra field, `mod_dir`, which
points to the location of the temporary working directory, so you should now
assign the result of calling `save_uwot` to the model you saved, e.g.
`model <- save_uwot(model, "my_model_file")`. This field is intended for use
with `unload_uwot`.
* `load_uwot` also returns the model with a `mod_dir` item for use with 
`unload_uwot`.
* `save_uwot` and `load_uwot` were not correctly handling relative paths.
* A previous bug fix to `load_uwot` in uwot 0.1.4 to work with newer versions
of RcppAnnoy (<https://github.com/jlmelville/uwot/issues/31>) failed in the 
typical case of a single metric for the nearest neighbor search using all 
available columns, giving an error message along the lines of:
`Error: index size <size> is not a multiple of vector size <size>`. This has now
been fixed, but required changes to both `save_uwot` and `load_uwot`, so 
existing saved models must be regenerated. Thank you to reporter 
[OuNao](https://github.com/OuNao).


# uwot 0.1.5

## Bug fixes and minor improvements

* The R API was being accessed from inside multi-threaded code to seed the 
(non-R) random number generators. Probably this was causing users in downstream
projects (seurat and monocle) to experience strange RcppParallel-related 
crashes. Thanks to [aldojongejan](https://github.com/aldojongejan)
for reporting this (<https://github.com/jlmelville/uwot/issues/39>).
* Passing a floating point value smaller than one to `n_threads` caused a crash.
This was particularly insidious if running with a system with only one default
thread available as the default `n_threads` becomes `0.5`. Now `n_threads` 
(and `n_sgd_threads`) are rounded to the nearest integer.
* Initialization of supervised UMAP should now be faster
(<https://github.com/jlmelville/uwot/issues/34>). Contributed by
[Aaron Lun](https://github.com/LTLA).

# uwot 0.1.4

## Bug fixes and minor improvements

* Fixed incorrect loading of Annoy indexes to be compatible with newer versions
of RcppAnnoy (<https://github.com/jlmelville/uwot/issues/31>). My thanks to
Dirk Eddelbuettel and Erik Bernhardsson for aid in identifying the problem.
* Fix for `ERROR: there is already an InterruptableProgressMonitor instance defined`.
* If `verbose = TRUE`, the `a`, `b` curve parameters are now logged.

# uwot 0.1.3

## Bug fixes and minor improvements

* Fixed an issue where the session would crash if the Annoy nearest neighbor
search was unable to find k neighbors for an item.

## Known issue

Even with a fix for the bug mentioned above, if the nearest neighbor index file
is larger than 2GB in size, Annoy may not be able to read the data back in. This
should only occur with very large or high-dimensional datasets. The nearest
neighbor search will fail under these conditions. A work-around is to set
`n_threads = 0`, because the index will not be written to disk and re-loaded
under these circumstances, at the cost of a longer search time. Alternatively,
set the `pca` parameter to reduce the dimensionality or lower `n_trees`, both of
which will reduce the size of the index on disk. However, either may lower the
accuracy of the nearest neighbor results.

# uwot 0.1.2

Initial CRAN release.

## New features

* New parameter, `tmpdir`, which allows the user to specify the temporary
directory where nearest neighbor indexes will be written during Annoy
nearest neighbor search. The default is `base::tempdir()`. Only used if
`n_threads > 1` and `nn_method = "annoy"`.

## Bug fixes and minor improvements

* Fixed an issue with `lvish` where there was an off-by-one error when
calculating input probabilities.

* Added a safe-guard to `lvish` to prevent the gaussian precision, beta,
becoming overly large when the binary search fails during perplexity
calibration.

* The `lvish` perplexity calibration uses the log-sum-exp trick to avoid
numeric underflow if beta becomes large.

# uwot 0.0.0.9010 (31 March 2019)

## New features

* New parameter: `pcg_rand`. If `TRUE` (the default), then a random number
generator from [the PCG family](https://www.pcg-random.org/) is used during the
stochastic optimization phase. The old PRNG, a direct translation of
an implementation of the Tausworthe "taus88" PRNG used in the Python
version of UMAP, can be obtained by setting `pcg_rand = FALSE`. The new PRNG is
slower, but is likely superior in its statistical randomness. This change in
behavior will be break backwards compatibility: you will now get slightly
different results even with the same seed.
* New parameter: `fast_sgd`. If `TRUE`, then the following combination of
parameters are set: `n_sgd_threads = "auto"`, `pcg_rand = FALSE` and `approx_pow
= TRUE`. These will result in a substantially faster optimization phase, at the
cost of being slightly less accurate and results not being exactly repeatable.
`fast_sgd = FALSE` by default but if you are only interested in visualization,
then `fast_sgd` gives perfectly good results. For more generic dimensionality
reduction and reproducibility, keep `fast_sgd = FALSE`.
* New parameter: `init_sdev` which specifies how large the standard deviation
of each column of the initial coordinates should be. This will scale any input
coordinates (including user-provided matrix coordinates). `init = "spca"` can
now be thought of as an alias of `init = "pca", init_sdev = 1e-4`. This may be
too aggressive scaling for some datasets. The typical UMAP spectral
initializations tend to result in standard deviations of around `2` to `5`, so
this might be more appropriate in some cases. If spectral initialization detects
multiple components in the affinity graph and falls back to scaled PCA, it
uses `init_sdev = 1`.
* As a result of adding `init_sdev`, the `init` options `sspectral`,
`slaplacian` and `snormlaplacian` have been removed (they weren't around for
very long anyway). You can get the same behavior by e.g.
`init = "spectral", init_sdev = 1e-4`. `init = "spca"` is sticking around
because I use it a lot.

## Bug fixes and minor improvements

* Spectral initialization (the default) was sometimes generating coordinates
that had too large a range, due to an erroneous scale factor that failed to
account for negative coordinate values. This could give rise to embeddings with
very noticeable outliers distant from the main clusters.
* Also during spectral initialization, the amount of noise being added had a
standard deviation an order of magnitude too large compared to the Python
implementation (this probably didn't make any difference though).
* If requesting a spectral initialization, but multiple disconnected components
are present, fall back to `init = "spca"`.
* Removed dependency on C++ `<random>` header. This breaks backwards
compatibility even if you set `pcg_rand = FALSE`.
* `metric = "cosine"` results were incorrectly using the unmodified Annoy
angular distance.
* Numeric matrix columns can be specified as the target for the `categorical`
metric (fixes <https://github.com/jlmelville/uwot/issues/20>).

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
* Internal changes and fixes thanks to a code review by
[Aaron Lun](https://github.com/ltla).

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
<https://github.com/jlmelville/uwot#mixed-data-types>.
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

* Embedding to `n_components = 1` was broken
(<https://github.com/jlmelville/uwot/issues/6>)
* User-supplied matrices to `init` parameter were being modified, in defiance of
basic R pass-by-copy semantics.

# uwot 0.0.0.9002 (August 14 2018)

## Bug fixes and minor improvements

* `metric = "cosine"` is working again for `n_threads` greater than `0`
(<https://github.com/jlmelville/uwot/issues/5>)

# uwot 0.0.0.9001

## New features

* *August 5 2018*. You can now use an existing embedding to add new points via
`umap_transform`. See the example section below.

* *August 1 2018*. Numerical vectors are now supported for supervised dimension
reduction.

* *July 31 2018*. (Very) initial support for supervised dimension reduction:
categorical data only at the moment. Pass in a factor vector (use `NA` for
unknown labels) as the `y` parameter and edges with bad (or unknown) labels are
down-weighted, hopefully leading to better separation of classes. This works
remarkably well for the Fashion MNIST dataset.

* *July 22 2018*. You can now use the cosine and Manhattan distances with the
Annoy nearest neighbor search, via `metric = "cosine"` and `metric =
"manhattan"`, respectively. Hamming distance is not supported because RcppAnnoy
doesn't yet support it.
