# Update from 0.7.1 to 0.7.2

fixed bug in `layout_with_stress3D` which only produced a 2D layout. No user facing changes.

## Downstream dependencies

All dependencies passed R CMD check

## Test environments
* ubuntu 20.04, R 4.1.2
* win-builder (devel and release)

## R CMD check results

There were no ERRORs, WARNINGs 

There was 1 NOTE on Ubuntu:

checking installed package size ... NOTE  
installed size is  7.1Mb  
sub-directories of 1Mb or more:  
      help   2.4Mb
      libs   4.5Mb
