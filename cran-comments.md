## New minor release

This is a new minor release to support a change in a reverse dependency which
has been asked to resubmit their package. To explain the situation: version
0.1.12 of uwot failed a reverse dependency check because of a broken vignette in
the bbknnR package. This was due to that package using a non-exported function
from this package, the internal structure of which changed when 0.1.12 was
created. I reversed that change as part of a successful re-submission creating
CRAN version 0.1.13.

This newest submission (0.1.14) exports a public function which bbknnR can use
in place of the non-exported function. I am aware of the policy of not
submitting new versions of packages to CRAN with such frequency. I apologize for
doing so, but request it be allowed to help the bbknnR package come back into
compliance with CRAN policies without undue delay.

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

There were two NOTEs:

N  checking installed package size ...
     installed size is 17.4Mb
     sub-directories of 1Mb or more:
       libs  16.9Mb

This is due to the underlying C++ implementation using templates.

The second NOTE:

* checking CRAN incoming feasibility ... [21s] NOTE
Maintainer: 'James Melville <jlmelville@gmail.com>'

Days since last update: 5

As mentioned above, this release is to provide a public function for the CRAN
package bbknnR to use, to replace its use of a non-exported internal function
from this package. Once this package is on CRAN, the bbknnR maintainer will be
able to submit a new version of their package and be back to following CRAN
policies. I regret causing extra work for the CRAN team but I think it's better
to give the bbknnR maintainer the opportunity to fix their package as soon as
possible. Uwe Ligges is aware of the issue with bbknnR and has asked the bbknnR
maintainers to resubmit their package, so I hope that this pattern of swift
resubmissions will be forgiven on this occasion.

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
