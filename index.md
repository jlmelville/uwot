# uwot

An R implementation of the [Uniform Manifold Approximation and
Projection (UMAP)](https://arxiv.org/abs/1802.03426) method for
dimensionality reduction of McInnes et al. (2018). Also included are the
supervised and metric (out-of-sample) learning extensions to the basic
method. Translated from the [Python
implementation](https://github.com/lmcinnes/umap).

## News

*August 29 2026* `uwot` version 0.2.5 has been released to CRAN. This is
to apply fixes for some minor issues surfaced by CRAN checks (transitive
headers in the C++ code caught by clang, a deprecation in test code) but
there are some other minor changes for nearest neighbor calculations.
See the
[NEWS.md](https://github.com/jlmelville/uwot/blob/master/NEWS.md) for
details.

## Installing

### From CRAN

``` r

install.packages("uwot")
```

### From GitHub

Building the development version from source requires a C++17 compiler.

- On Windows, install the [Rtools version for your R
  release](https://cran.r-project.org/bin/windows/Rtools/).
- On macOS, install the Xcode Command Line Tools with
  `xcode-select --install`. The [R for macOS tools
  page](https://mac.r-project.org/tools/) lists compilers for
  dependencies built from source.
- On Linux, install R’s development headers, `make`, and a C++17
  compiler through your distribution.

``` r

install.packages("pak")
pak::pak("jlmelville/uwot")
```

## Quick start

Embed the four measurements in R’s built-in `iris` data, then plot one
point per flower. The plotting symbols identify the species.

``` r

library(uwot)

iris_umap <- umap2(iris[, 1:4], seed = 42)
species <- iris$Species
plot(
  iris_umap,
  pch = as.integer(species),
  xlab = "UMAP 1",
  ylab = "UMAP 2"
)
legend(
  "topright",
  legend = levels(species),
  pch = seq_along(levels(species)),
  bty = "n"
)
```

## MNIST example

To embed the 70,000 handwritten digits in MNIST, first install the
`snedata` package from GitHub.

``` r

install.packages("pak")
pak::pak("jlmelville/snedata")
```

Download the data, fit UMAP to the pixel values, and color each point by
its digit label.

``` r

mnist <- snedata::download_mnist()
mnist_umap <- uwot::umap(mnist, min_dist = 0.001, verbose = TRUE)
digit_colors <- grDevices::adjustcolor(rainbow(nlevels(mnist$Label)), alpha.f = 0.1)
plot(
  mnist_umap,
  cex = 0.1,
  col = digit_colors[mnist$Label],
  main = "R uwot::umap",
  xlab = "",
  ylab = ""
)
```

![UMAP embedding of MNIST, colored by
digit](reference/figures/mnist-r.png)

UMAP embedding of MNIST, colored by digit

## Optional neighbor backends

For approximate nearest-neighbor search, install the package for the
method you want to use.

``` r

install.packages(c("RcppHNSW", "rnndescent"))
```

Select HNSW or NN-Descent with `nn_method`. These examples use the
built-in iris measurements.

``` r

iris_hnsw <- uwot::umap(iris[, 1:4], nn_method = "hnsw")
iris_nnd <- uwot::umap(iris[, 1:4], nn_method = "nndescent")
```

## Documentation

<https://jlmelville.github.io/uwot/>. For more examples see the [get
started](https://jlmelville.github.io/uwot/articles/uwot.html) doc.
There are plenty of
[articles](https://jlmelville.github.io/uwot/articles/index.html)
describing various aspects of the package.

## License

[GPLv3 or later](https://www.gnu.org/licenses/gpl-3.0.txt).

## Citation

If you want to cite the use of uwot, then use the output of running
`citation("uwot")` (you can do this with any R package).

## See Also

- The [UMAP reference implementation](https://github.com/lmcinnes/umap)
  and [publication](https://arxiv.org/abs/1802.03426).
- The [UMAP R package](https://cran.r-project.org/package=umap) (see
  also its [github repo](https://github.com/tkonopka/umap)), predates
  `uwot`’s arrival on CRAN.
- Another R package is
  [umapr](https://github.com/ropensci-archive/umapr), but it is no
  longer being maintained.
- [umappp](https://github.com/libscran/umappp) is a full C++
  implementation, and [yaumap](https://github.com/LTLA/yaumap) provides
  an R wrapper. The batch implementation in umappp are the basis for
  uwot’s attempt at the same.
- `uwot` uses the
  [RcppProgress](https://cran.r-project.org/package=RcppProgress)
  package to show a text-based progress bar when `verbose = TRUE`.
