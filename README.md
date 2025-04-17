
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ncappc: NCA Calculations and population model diagnosis

<!-- badges: start -->

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ncappc)](https://CRAN.R-project.org/package=ncappc)
[![R-CMD-check](https://github.com/UUPharmacometrics/ncappc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/UUPharmacometrics/ncappc/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/UUPharmacometrics/ncappc/graph/badge.svg)](https://app.codecov.io/gh/UUPharmacometrics/ncappc)
<!-- badges: end -->

ncappc performs NCA Calculations and population model diagnosis using
posterior predictive checks generated from data simulated by a
population model.

ncappc is a flexible tool that can perform:

1.  Traditional non-compartmental analysis (NCA) and
2.  Simulation-based posterior predictive checks for population
    pharmacokinetic (PK) and/or pharmacodynamic (PKPD) models using NCA
    metrics.

You can read more at the [website for the stable version of
ncappc](https://uupharmacometrics.github.io/ncappc/) or the [website for
the development version of
ncappc](https://uupharmacometrics.github.io/ncappc/dev/)

## Installation

You need to have R installed. Download the latest version of R from
www.r-project.org.

You can install the latest stable release from CRAN:

``` r
install.packages("ncappc")
```

You can install the development version of ncappc from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("UUPharmacometrics/ncappc")
```

<!-- ## Example -->

<!-- This is a basic example which shows you how to solve a common problem: -->

<!-- ```{r example} -->

<!-- library(ncappc) -->

<!-- ## basic example code -->

<!-- ``` -->

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->

<!-- summary(cars) -->

<!-- ``` -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. -->

<!-- You can also embed plots, for example: -->

<!-- ```{r pressure, echo = FALSE} -->

<!-- plot(pressure) -->

<!-- ``` -->

<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
