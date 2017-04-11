
<!-- README.md is generated from README.Rmd. Please edit that file -->
**NOTE: This package is still under very active development**

[![Travis-CI Build Status](https://travis-ci.org/stephenslab/rssr.svg?branch=master)](https://travis-ci.org/stephenslab/rssr)

Installation
------------

The first step in installation is to ensure that the `devtools` library is installed

``` r
install.packages('devtools')
devtools::install_github("stephenslab/rssr",ref="v0.1.1-alpha",build_vignettes = TRUE)
```

Usage
-----

The functions `rss_varbvsr_naive` and `rss_varbvsr_squarem`contain the default and accelerated implementations of RSS using variational bayes. A vignette can be found after installation by running `vignette("example")`
