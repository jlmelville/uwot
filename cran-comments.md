## Patch release

This is a patch release to fix an ERROR that has occurred on the CRAN checks
with the previous release.

## Recent Resubmission

The previous version was submitted to CRAN less than a week ago. I regret the
inconvenience of having to resubmit so soon, but on balance decided that fixing
an ERROR status in CRAN checks as soon as possible was the best course of
action.

## Test environments

* local ubuntu 23.10 R 4.3.1
* ubuntu 22.04 (on github actions), R 4.2.3, R 4.3.3, devel
* ubuntu 22.04 (on rhub) devel clang-ASAN
* ubuntu 22.04 (on rhub) R 4.4.0 beta gcc
* Fedora 38 (on rhub) devel valgrind
* win-builder (devel)
* local Windows 11 build, R 4.3.3
* Windows Server 2022 (on github actions), R 4.2.3, R 4.3.3
* Windows Server 2012 (on appveyor) R 4.3.3
* local mac OS X Sonoma R 4.3.3
* mac OS X Monterey (on github actions) R 4.3.3

## R CMD check results

There are two ERRORs on r-release-macos-arm64: 

Check: tests
Check: re-building of vignette outputs

In both cases the source of the error is the same and is due to an ABI change in
a version of the `Matrix` package interacting poorly with the `irlba` package. I
have made changes to avoid the erroneous code path being encountered. This
submission is intended to fix these errors.

There are no WARNINGs.

There is 1 NOTE remaining from the previous release on r-prerel-macos-arm64, 
r-prerel-macos-x86_64, r-release-macos-arm64, r-release-macos-x86_64, 
r-oldrel-macos-arm64:

Check: installed package size
Result: NOTE 
    installed size is 13.8Mb
    sub-directories of 1Mb or more:
      libs  12.8Mb

This is due to the underlying C++ implementation using templates.

## revdepcheck results

We checked 65 reverse dependencies (34 from CRAN + 31 from Bioconductor), 
comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
