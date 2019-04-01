## Release Summary

* This is a new release.

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

There was one NOTE:

* checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.

This is expected due to linking to the RcppParallel package.

There was a message about possibly mis-spelled wordss in DESCRIPTION:
  
     Healy (7:16)
     McInnes (7:4)
     SNE (8:75)
     UMAP (2:59, 6:50)
     arxiv (7:30)
     
These are spelled correctly.

## Downstream dependencies

None.
