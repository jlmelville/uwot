## Release Summary

This is a patch release to maintain compatibility with a forthcoming version of
the RcppAnnoy dependency.

## Test environments

* ubuntu 14.04 (on travis-ci), R 3.4.4, R 3.6.0, R-devel
* ubuntu 16.04 (on rhub), R 3.6.1
* fedora 30 (on rhub), R-devel
* mac OS X High Sierra (on travis-ci), R 3.5.3, R 3.6.1
* local Windows 10 build, R 3.6.1
* Windows Server 2008 (on rhub) R-devel
* Windows Server 2012 (on appveyor) R 3.6.1
* win-builder (devel)

## R CMD check results

There were no ERRORs or WARNINGs.

There was one NOTE:

* checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.

This is expected due to linking to the RcppParallel package.

With r-hub checking on Windows only there was also:

"N  checking for non-standard things in the check directory
   Found the following files/directories:
     'examples_x64' 'tests_i386' 'tests_x64'
     'RcppHNSW-Ex_i386.Rout' 'RcppHNSW-Ex_x64.Rout' 'examples_i386'"

This would seem to be something to do with r-hub rather than a real problem.

There was a message about possibly mis-spelled words in DESCRIPTION:
  
  McInnes (7:42)
  Rcpp (11:73)
  UMAP (2:59)
  al (7:53, 10:38)
  et (7:50, 10:35)
  uwot (12:57)
     
These are spelled correctly.

## CRAN checks

A gcc-UBSAN issue is reported, due to a library used by the RcppParallel 
package. The RcppParallel package CRAN check also reports this issue, so I don't
think it's something that I can fix in this package.

## Downstream dependencies

There are 4 downstream dependencies. 

* 'dyndimred' fails R CMD CHECK because of test failures when it attempts to 
install the 'destiny' package, which is not available on R 3.6.1. This is 
unrelated to the new version of uwot. When these tests are disabled, the other
tests, and R CMD CHECK as a whole, completes successfully.

* The other three dependencies completed R CMD CHECK without issues.
