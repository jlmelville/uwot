## Resubmission

This is a resubmission. In this version I have:

* Resolved usages of dependencies that resulted in UBSAN errors.

## Test environments

* ubuntu 16.04 (on travis-ci), R 3.5.3, R 3.6.2, R-devel
* ubuntu 16.04 (on rhub), R 3.6.1
* fedora 31 (on rhub), R-devel
* mac OS X High Sierra (on travis-ci), R 3.5.3, R 3.6.2
* local Windows 10 build, R 3.6.2, R 3.6.3
* local Pop OS (Linux) build R 3.6.3
* Windows Server 2008 (on rhub) R-devel
* Windows Server 2012 (on appveyor) R 3.6.2

## R CMD check results

There were no ERRORs or WARNINGs.

* For r-hub with Windows Server 2008 R2 SP1, R-devel, 32/64 bit there was one 
  NOTE:

checking for non-standard things in the check directory ... NOTE
  'examples_i386' 'examples_x64' 'tests_i386' 'tests_x64'
Found the following files/directories:
  'uwot-Ex_i386.Rout' 'uwot-Ex_x64.Rout'

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

The following issues were reported in the previous submission:

* r-patched-solaris-x86 had an ERROR, because the RcppAnnoy package was 
unavailable.

* A gcc-UBSAN issue was reported due to libraries used by RcppAnnoy and 
RcppParallel.

Both should now be fixed.

## Downstream dependencies

There are 4 downstream CRAN dependencies and 3 Bioconductor dependencies:

* Five of these packages report either NOTEs or a mixture of WARNINGs and NOTEs.
None are related to 'uwot'.
* The other two dependencies completed R CMD CHECK without issues.
