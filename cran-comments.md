## Changes

Changes in this version of ncappc are:
  

## Test environments
* local macOS (10.13.6) install, R release
* ubuntu 14.04 (on travis-ci), R release and devel version
* Windows Server 2012 R2 x64 (on AppVeyor), R release 
* win-builder, R release and devel version

## R CMD check results
For macOS, travis and AppVeyor there were no ERRORs, WARNINGs or NOTEs. 

For win-builder release and devel versions there was one note:
  
  * checking CRAN incoming feasibility ... NOTE
+ Maintainer: 'Andrew C. Hooker <andrew.hooker@farmbio.uu.se>'
+ Possibly mis-spelled words in DESCRIPTION:
  pharmacometric (19:9)

The first portion is just stating that I am the maintainer.  

The second portion does not recognize the well known word "pharmacometric",
see for example: https://en.wikipedia.org/wiki/Pharmacometrics

## Downstream dependencies
There are currently no downstream dependencies for this package.
