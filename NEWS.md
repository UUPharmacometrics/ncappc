# ncappc (development version)

* This is a major release (to version 1.0.0), since the package 
  is in a stable state and is not considered experimental anymore. 
  
* Created a webpage for the package using pkgdown. See
  https://uupharmacometrics.github.io/ncappc/
  
* Updated which time and concentration values are used for 
  extrapolation for lamda_z when interpolation occurs (#5).
  
* Updated tidyverse functions that have been depreciated. 

* Removed the Cairo package from the dependencies of this package.

# ncappc 0.3.0

* New mean and variance of NCA metrics histogram plots

* Updating how reports are generated and standardizing how they appear

* Change how estimate_nca works, making tidy calculations the default
  and assuring that that the same function works for both real and simulated
  data

* Added parallel computation of NCA metrics for simulated data

* Added a number of automatic tests used for checking the package 
  including results from ncappc publication.

* Added a `NEWS.md` file to track changes to the package.

* Fixed depencies between ggplot2, grid and scales

* Updated vignettes

* More examples added to documentation

* Allowing extrapolation for fewer than three points in the elimination phase
  of a profile

* Allow timing messages to be turned on and off.

* Bug fixes

* linear-log option was changed to linearup-logdown 

* Central tendency measure for NCA metrics changed from mean to median

* Updated documentation

* Changed package maintainer

* MRTINF pitfall fixed

* Added continuous integration

* Made onlyNCA work as it should




