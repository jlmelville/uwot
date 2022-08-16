## Resubmission

This is a re-submission to deal with a broken vignette in a CRAN reverse
dependency which is using a non-public function.

The original submission was a new feature and bug fix release to maintain 
compatibility with a forthcoming version of the Matrix package.

## Test environments

* mac OS X Big Sur (on github actions), R 4.2.1
* Fedora Linux, R-devel, clang, gfortran (on R-hub)
* Ubuntu Linux 20.04.1 LTS, R-release, GCC (on R-hub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (on R-hub)
* ubuntu 22.04 R 4.2.1 (with valgrind)
* ubuntu 20.04 (on github actions) R 4.1.3, R 4.2.1, devel
* Windows Server 2022, R-devel, 64 bit (on R-hub)
* win-builder R devel
* local Windows 11 build, R 4.2.1
* Windows Server 2012 (on appveyor) R 4.2.1
* Windows Server 2022 (on github actions) R 4.1.3, 4.2.1

## R CMD check results

There are no WARNINGs or ERRORs.

There was one NOTE:

N  checking installed package size ...
     installed size is 17.4Mb
     sub-directories of 1Mb or more:
       libs  16.9Mb

This is due to the underlying C++ implementation using templates.

* There was a message about possibly mis-spelled words in DESCRIPTION:
  
  McInnes (7:42)
  Rcpp (11:73)
  UMAP (2:59)
  al (7:53, 10:38)
  et (7:50, 10:35)
  uwot (12:57)
     
These are spelled correctly.

## Downstream dependencies

We checked 54 reverse dependencies (28 from CRAN + 26 from Bioconductor), 
comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
