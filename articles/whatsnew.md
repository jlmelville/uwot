# What's New

A place for me to put the old less-structured updates that I post at the
top of the `README.md` file for the package. When they aren’t new any
more, they will get moved here. You should look at the
[Changelog](https://jlmelville.github.io/uwot/news/index.html) for
fuller details.

*November 10 2025* `uwot` version 0.2.4 has been released to CRAN. This
was mainly to avoid a test potentially starting to fail due to an
incorrect use of `testthat`. However, as part of some other small fixes,
optional dependencies (e.g. RSpectra) will now be correctly detected and
used if installed even if not loaded. This could have an effect on
output (but the old behavior was a bug).

*February 24 2025* `uwot` version 0.2.3 has been released to CRAN. This
release mainly fixes some bugs, including one that was causing an error
with an upcoming version of R-devel. One new feature: set
`rng_type = "deterministic"` to use a deterministic sampling of vertices
during the optimization phase which will give faster and more
reproducible output at the cost of accuracy. The idea for this came
straight from [Leland McInnes via
Reddit](https://www.reddit.com/r/MachineLearning/comments/1gsjfq9/comment/lxip9wy/).

*April 21 2024* As ordained by prophecy, version 0.2.2 of `uwot` has
been released to CRAN. `RSpectra` is back as a main dependency and I
*thought* I had worked out a clever scheme to avoid the failures seen in
some installations with the `irlba`/`Matrix` interactions. This releases
fixes the problem on all the systems I have access to (including GitHub
Actions CI) but some CRAN checks remain failing. How embarrassing. That
said, if you have had issues, it’s possible this new release will help
you too.

*April 18 2024* Version 0.2.1 of `uwot` has been released to CRAN. Some
features to be aware of:
[RcppHNSW](https://cran.r-project.org/package=rnndescent) and
[rnndescent](https://cran.r-project.org/package=rnndescent) are now
supported as optional dependencies. If you install and load them, you
can use them as an alternative to RcppAnnoy in the nearest neighbor
search and should be faster. Also, a new `umap2` function has been
added, with updated defaults compared to `umap`. Please see the updated
and new articles on
[HNSW](https://jlmelville.github.io/uwot/articles/hnsw-umap.html),
[rnndescent](https://jlmelville.github.io/uwot/articles/rnndescent-umap.html),
[working with sparse
data](https://jlmelville.github.io/uwot/articles/sparse-data-example.html)
and [umap2](https://jlmelville.github.io/uwot/articles/umap2.html). I
consider this worthy of moving from `0.1.x` to `0.2.x`, but in the
interests of full disclosure, on-going [irlba
problems](https://github.com/jlmelville/uwot/issues/115) has caused a
CRAN check failure, so we might be onto 0.2.2 sooner than I’d like.

*November 26 2023* Happy 1 million CRAN downloads to `uwot`. To
celebrate (actually it’s just a coincidence) I have reorganized the
horror-show that was the sprawling README into actual articles, via the
fantastic [pkgdown](https://pkgdown.r-lib.org/) package. Find it at
<https://jlmelville.github.io/uwot/>.

*June 28 2023* Version 0.1.16 has been released to CRAN. This is a very
minor tweak to 0.1.15 to further support the new release of
[RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy).

*June 28 2023* Version 0.1.16 has been released to CRAN. This is a very
minor tweak to 0.1.15 to further support the new release of
[RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy).

*June 26 2023* Version 0.1.15 has been released to CRAN. This is to
support a new release of
[RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy), but there are
also some bug fixes and other minor improvements. There are some new
functions: `optimize_graph_layout` will carry out the UMAP optimization
step on a sparse similiarity matrix, e.g. the output of
`similarity_graph`. `simplicial_set_union` and
`simplicial_set_intersect` provide ways to merge different views of the
same data into one sparse similiarity matrix. As usual, `NEWS.md` has
all the details.

*August 22 2022* Just when you least expected it, version 0.1.14 has
been released to CRAN (the `NEWS` file on CRAN calls it `0.1.13.9000`
because I forgot to update that file, but let’s keep that amongst
ourselves). This release includes a bug fix for `umap_transform` when
you use external nearest neighbors and new function `similarity_graph`,
to support extracting just the high dimensional fuzzy simplicial set.

*August 16 2022* Version 0.1.13 has been released to CRAN (0.1.12 was a
failed submission). Among other things you can now pass your own nearest
neighbors data in sparse matrix form. Also there is an option to
reproduce relative cluster density by [approximating the densMAP
method](https://jlmelville.github.io/uwot/articles/leopold.html). See
the
[NEWS](https://github.com/jlmelville/uwot/blob/master/NEWS.md#uwot-0113)
page for more.

*December 12 2021* Version 0.1.11 has been released to CRAN. It is now
possible to get reproducible results (for a given value of `set.seed`)
when running the optimization step with multiple threads
(`n_sgd_threads` greater than 1). You may need to increase `n_epochs` to
get similar levels of convergence. To run in this mode, set
`batch = TRUE`. Thanks to [Aaron Lun](https://github.com/LTLA) who came
up with the design for this and also implemented it in his [umappp C++
library](https://github.com/libscran/umappp). See `NEWS.md` for other
changes.

*December 15 2020* Version 0.1.10 has been released to CRAN. This is
mainly to maintain compatibility with RcppAnnoy, but also a small change
was made to avoid it grinding away pointlessly in the presence of `NA`
values, based on an observation by David McGaughey on Twitter (which I
can no longer link to).

*November 15 2020* Version 0.1.9 has been released to CRAN. The main
addition is support for the Pearson correlation. Also, a slight license
change from GPL-3 to GPL-3 or later.

*August 1 2020* New metric supported: Pearson correlation (with
`metric = "correlation"`). This should give similar results to the
Python UMAP (and sklearn) implementation of the `correlation` metric.

*March 16 2020* A new version (0.1.8) is on CRAN. This is a minor
release in terms of features, but you can now export the UMAP graph
(<https://github.com/jlmelville/uwot/issues/47>), and there are some bug
fixes for: loading Annoy indexes
(<https://github.com/jlmelville/uwot/issues/31>), reproducibility across
platforms (<https://github.com/jlmelville/uwot/issues/46>) and we no
longer use RcppParallel for the multi-threading support, which should
lead to fewer installation problems.

*March 4 2020* I had to cancel my submission of version 0.1.7 to CRAN
because of a broken example in a library using uwot. In the mean time I
have switched to using `std::thread` rather than tinythread++.

*March 1 2020* Version 0.1.6 was rejected from CRAN due to undefined
behavior issues that originate from RcppAnnoy and RcppParallel. I am
hopeful that the Annoy behavior is fixed and a suitable version of
RcppAnnoy will be released onto CRAN eventually. The RcppParallel issues
originate with the use of [tbb](https://github.com/oneapi-src/oneTBB)
and seems much harder to deal with. As there is no way to use
RcppParallel without tbb yet, I am temporarily replacing the use of
RcppParallel with just a subset of the code needed to run parallel for
loops with the [tinythread++](https://tinythreadpp.bitsnbites.eu/)
library.

*December 4 2019* Version 0.1.5 released on CRAN. This fixes a couple of
crash bugs, including one where the R API was being called from inside a
thread. This may have been causing the issues seen by users of
[monocle](https://github.com/cole-trapnell-lab/monocle3/issues/186) and
[seurat](https://github.com/satijalab/seurat/issues/2256).

*September 23 2019* Version 0.1.4 released on CRAN. This ensures
compatibility with RcppAnnoy 0.0.13 when using `load_uwot`.

*April 6 2019*. uwot is now on
[CRAN](https://cran.r-project.org/package=uwot). Also, some
minor-to-horrible bugs in the `lvish` perplexity routine have been
fixed.

For visualization purposes, it seems reasonable to use the old PRNG
(`pcg_rand = FALSE`), along with multiple threads during SGD
(`n_sgd_threads = "auto"`), and the UMAP gradient approximation
(`approx_pow = TRUE`), which combined will show a very noticeable speed
up during optimization. I have added a new parameter, `fast_sgd`, which
if set to `TRUE`, sets these options for you.
