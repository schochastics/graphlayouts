# Update from 0.2.0 to 0.3.0

This submission contains new layout functions (pivot MDS, sparse stress and dynamic layouts)
and a series of small fixes to make the package integrate well with `ggraph`.

# Test environments
* ubuntu 18.04, R 3.6.1
* win-builder (devel and release)
* macOS 10.11 El Capitan, R-release

## R CMD check results

There were no ERRORs, WARNINGs 

There was 1 NOTE on Ubuntu:

checking installed package size ... NOTE  
  installed size is  5.2Mb  
  sub-directories of 1Mb or more:  
    help   1.9Mb  
    libs   3.2Mb  

RcppArmadillo was added in this version for two layout algorithms which appeared to
increase the installation size.

## Downstream dependencies

The only package that depends on `graphlayouts` is the `SNAhelper`, which passed 
R CMD check

# Update from 0.1.0 to 0.2.0 (Resubmission 1)

The `README.md` contained an invalid URL (CRAN URL not in canonical form)
This is fixed now

# Update from 0.1.0 to 0.2.0

This update includes some bugfixes and a faster implementation for `layout_as_backbone()`

## Test environments

* ubuntu 18.04, R 3.6.0
* ubuntu 16.04 (on travis-ci), R 3.6.0
* win-builder (devel and release)
* macOS 10.11 El Capitan, R-release

## R CMD check results

There were no ERRORs, WARNINGs, or NOTEs 

## Downstream dependencies

The only package that depends on `graphlayouts` is the `SNAhelper`, which passed 
R CMD check


# Resubmission of inital submission

## Test environments

* ubuntu 18.04, R 3.5.2
* ubuntu 14.04 (on travis-ci), R 3.5.2
* win-builder (devel and release)
* macOS 10.11 El Capitan, R-release

## R CMD check results

There were no ERRORs, WARNINGs 

There was 1 NOTE:

* This is a new release.

## Comments

Thanks for your feedback on the resubmission. All examples were unwrapped and
no dontrun or donttest is included anymore. examples run in 3.9s.

