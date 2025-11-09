## Patch release

This is a patch release for compatibility with a new release of the `testthat`
package, plus some other small bug fixes.

## Test environments

* local ubuntu 25.10 R 4.5.1
* ubuntu 24.04 (on github actions), R 4.4.3, R 4.5.2, devel
* ubuntu 22.04 (on rhub) devel clang-ASAN
* ubuntu 22.04 (on rhub) devel gcc13
* Fedora 38 (on rhub) devel valgrind
* win-builder (devel)
* local Windows 11 build, R 4.5.2
* Windows Server 2022 (on github actions), R 4.4.3, R 4.5.2
* Windows Server 2012 (on appveyor) R 4.5.2
* local mac OS X Sequoia R 4.5.2
* mac OS X Sequoia (on github actions) R 4.5.2

## R CMD check results

There is 1 NOTE:

Version: 0.2.3
Check: installed package size
Result: NOTE 
    installed size is 16.6Mb
    sub-directories of 1Mb or more:
      libs  15.6Mb
Flavors: r-oldrel-macos-arm64, r-oldrel-macos-x86_64

This is due to the underlying C++ implementation using templates.

## revdepcheck results

We checked 76 reverse dependencies (39 from CRAN + 37 from Bioconductor), 
comparing R CMD check results across CRAN and dev versions of this package.

* We saw 0 new problems
* We failed to check 0 packages
