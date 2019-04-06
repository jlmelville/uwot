## Resubmission

This is a resubmission. In this version I have:

* Replaced double-blanks in Description field of the DESCRIPTION with four
spaces.

* Added single quotes to 'Rcpp' in the Description field of the DESCRIPTION.

* Added a URL to the uwot website in the Description field of the DESCRIPTION.

* Replaced \dontrun with \donttest in examples in the Rd files.

* Added a tmpdir parameter to functions to explicitly secify the temporary 
directory used to store temporary data.

## Resubmission

This is a resubmission. In this version I have:

* Expanded the Description field of the DESCRIPTION file, including explaining
the UMAP acronym and corrected the arXiv reference.

* Provided runnable examples in the .Rd files for every exported function.

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
