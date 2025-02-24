## Patch release

This is a patch release to fix an ERROR that occurs with an upcoming change in
R-devel.

## Test environments

* local ubuntu 24.10 R 4.4.1
* ubuntu 24.04 (on github actions), R 4.3.3, R 4.4.2, devel
* ubuntu 22.04 (on rhub) devel clang-ASAN
* ubuntu 22.04 (on rhub) R 4.4.3 beta gcc
* Fedora 38 (on rhub) devel valgrind
* win-builder (devel)
* local Windows 11 build, R 4.4.2
* Windows Server 2022 (on github actions), R 4.3.3, R 4.4.2
* Windows Server 2012 (on appveyor) R 4.4.2
* local mac OS X Sequoia R 4.3.3
* mac OS X Sonoma (on github actions) R 4.4.2

## R CMD check results

There is one ERROR on r-devel-linux-x86_64-debian-gcc, 
r-devel-linux-x86_64-fedora-clang and r-devel-linux-x86_64-fedora-gcc:

    ══ Failed tests ════════════════════════════════════════════════════════════════
    ── Failure ('test_errors.R:28:1'): (code run outside of `test_that()`) ─────────
    `umap(iris10, pca = 500)` threw an error with unexpected message.
    Expected match: "'pca' must be <="
    Actual message: "'length = 40' in coercion to 'logical(1)'"

This submission is intended to fix this error.

There is one ERROR on r-oldrel-macos-arm64:

Error in `irlba::irlba(L, nv = n, nu = 0, maxit = iters)`: function 'as_cholmod_sparse' not provided by package 'Matrix'

This is due to an ABI change in a version of the `Matrix` package interacting
poorly with the `irlba` package. I am no longer able to reproduce this error on
any of the platforms I have tested on. I believe that recompiling the `irlba`
(and maybe `Matrix`?) package from source would fix the issue. More details
are at <https://github.com/bwlewis/irlba/issues/70>. At any rate I do not 
believe this is a problem with my package.

There are no WARNINGs.

There is 1 NOTE:

Check: installed package size
Result: NOTE 
    installed size is 13.8Mb
    sub-directories of 1Mb or more:
      libs  12.8Mb
Flavors: r-release-macos-arm64, r-release-macos-x86_64, r-oldrel-macos-arm64, r-oldrel-macos-x86_64

This is due to the underlying C++ implementation using templates.

## revdepcheck results

We checked 69 reverse dependencies (36 from CRAN + 33 from Bioconductor), 
comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
