## New minor release

This is a new minor release to support a forthcoming release of one of this
package's dependencies, RcppAnnoy. Testing of the forthcoming submission of
RcppAnnoy with the previous version of this package identified the need for a
further modification of uwot to support one of uwot's reverse dependencies.
Reverse dependency checking of uwot could not have found this until RcppAnnoy is
updated on CRAN.

## Recent resubmission

0.1.15 was submitted only two days ago. My apologies for the swift resubmission.
This is to account for the issues found as noted above.

## Test environments

* mac OS X Monterey (on github actions), R 4.3.1
* Fedora Linux, R-devel, clang, gfortran (on R-hub)
* Debian, R-release, GCC (on R-hub)
* ubuntu 23.04 R 4.2.2 (with valgrind)
* ubuntu 22.04 (on github actions) R 4.2.3, R 4.3.1, devel
* Windows Server 2022, R-devel, 64 bit (on R-hub)
* win-builder R devel
* local Windows 11 build, R 4.3.1
* Windows Server 2012 (on appveyor) R 4.3.1
* Windows Server 2022 (on github actions) R 4.2.3, 4.3.1

## R CMD check results

There are no WARNINGs or ERRORs.

There is one NOTE:

* checking installed package size ... NOTE
  installed size is 18.0Mb
  sub-directories of 1Mb or more:
    libs  17.5Mb

This is due to the underlying C++ implementation using templates.

## Package Check Problems

These issues were fixed with the previous release of uwot, but as not all checks
are complete yet, I am including this here to avoid confusion:

* There is a NOTE: `Specified C++11: please drop specification unless essential`
-- this has been fixed with this submission.

* There is an `M1mac` test failure <https://www.stats.ox.ac.uk/pub/bdr/M1mac/uwot.out>
-- this has been fixed with this submission.

## Downstream dependencies

We checked all 32 CRAN reverse dependencies, comparing R CMD check results 
across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
