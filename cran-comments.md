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

