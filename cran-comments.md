## Resubmission

This is a resubmission to reformat arxiv links to doi format.

## Test environments

* local ubuntu 23.10 R 4.3.1
* ubuntu 22.04 (on github actions), R 4.3.3, devel
* ubuntu 22.04 (on rhub) devel clang-ASAN
* ubuntu 22.04 (on rhub) devel clang18
* ubuntu 22.04 (on rhub) R 4.4.0 beta gcc
* Fedora 38 (on rhub) devel valgrind
* win-builder (devel)
* local Windows 11 build, R 4.3.3
* Windows Server 2022 (on github actions), R 4.3.3
* Windows Server 2012 (on appveyor) R 4.3.3
* local mac OS X Sonoma R 4.3.3
* mac OS X Monterey (on github actions) R 4.3.3

## R CMD check results

There are no WARNINGs or ERRORs.

There is 1 NOTE remaining from the previous release on r-release-macos-arm64, 
r-release-macos-x86_64, r-oldrel-macos-arm64:

Check: installed package size
Result: NOTE 
    installed size is 18.7Mb
    sub-directories of 1Mb or more:
      libs  18.1Mb

This is due to the underlying C++ implementation using templates.

## revdepcheck results

We checked 65 reverse dependencies (34 from CRAN + 31 from Bioconductor), 
comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
