## Test environments

* ubuntu 16.04 (on travis-ci), R 3.5.3, R 3.6.2, R-devel
* ubuntu 16.04 (on rhub), R 3.6.1
* fedora 31 (on rhub), R-devel
* mac OS X High Sierra (on travis-ci), R 3.5.3, R 3.6.2
* local Windows 10 build, R 3.6.2
* Windows Server 2008 (on rhub) R-devel
* Windows Server 2012 (on appveyor) R 3.6.2

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

r-patched-solaris-x86	has an ERROR. This is because the RcppAnnoy package is
unavailable.

### UBSAN issue

A gcc-UBSAN issue is reported. This originates from the dependency RcppAnnoy and 
is ultimately due to how memory management is implemented in the Annoy 
library that RcppAnnoy wraps. This is a purposeful design decision in Annoy 
and not something that can be fixed in this package (or in RcppAnnoy).

## Downstream dependencies

There are 4 downstream CRAN dependencies and 3 Bioconductor dependencies:

* Five of these packages report either NOTEs or a mixture of WARNINGs and NOTEs.
None are related to 'uwot'.
* The other two dependencies completed R CMD CHECK without issues.
