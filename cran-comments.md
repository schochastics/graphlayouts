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

Thanks for your feedback on the initial submission. Below are my comments on your suggestions:

- Package names in Description field of DESCRIPTION were surounded by "'"
- A reference was added in the Description field of DESCRIPTION
- examples unwrapped from dontrun when feasible and keeping the runtime under 5s (layout_as_backbone example takes too long)
- more examples added 
- It was pointed out that I should add all authors and copyright holders in the Authors@R field with 
the appropriate roles. The complete package (including Rcpp code) was written by myself. 
If you suspect any c/p code or missing authorship of code please indicate the piece of code and I will look into it. 

