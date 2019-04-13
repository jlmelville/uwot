## CRAN checks

A gcc-UBSAN issue is reported, due to a library used by the RcppParallel 
package. The RcppParallel package CRAN check also reports this issue, so I don't
think it's something that I can fix in this package.

## Release Summary

This is a patch release for a bug in C++ nearest neighbor code that could cause
the session to terminate.

## Test environments

* ubuntu 14.04 (on travis-ci), R 3.4.4, R 3.5.2, R-devel
* ubuntu 16.04 (on rhub), R 3.4.4
* debian (on rhub), R-devel
* mac OS X High Sierra (on travis-ci), R 3.4.4, R 3.5.3
* local Windows 10 build, R 3.5.3
* Windows Server 2008 (on rhub) R-devel
* Windows Server 2012 (on appveyor) R 3.5.2
* win-builder (devel)

## R CMD check results

There were no ERRORs or WARNINGs.

There were two NOTEs:

* checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.

This is expected due to linking to the RcppParallel package.

* Days since last update: 1

I apologize for the submission of a new release so soon after initial
acceptance. The bug is unlikely to occur in practice but because it causes the
session to crash, I considered it urgent to fix it ASAP.

There was a message about possibly mis-spelled words in DESCRIPTION:
  
  McInnes (7:42)
  Rcpp (11:73)
  UMAP (2:59)
  al (7:53, 10:38)
  et (7:50, 10:35)
  uwot (12:57)
     
These are spelled correctly.

## Downstream dependencies

None.
