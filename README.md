# AtmoSwing R toolbox [![DOI](https://zenodo.org/badge/90713710.svg)](https://zenodo.org/badge/latestdoi/90713710) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/atmoswing/R-toolbox?branch=master&svg=true)](https://ci.appveyor.com/project/atmoswing/R-toolbox) [![Build Status](https://travis-ci.org/atmoswing/r-toolbox.svg?branch=master)](https://travis-ci.org/atmoswing/r-toolbox) [![Coverage Status](https://coveralls.io/repos/github/atmoswing/R-toolbox/badge.svg?branch=master)](https://coveralls.io/github/atmoswing/R-toolbox?branch=master)
Postprocessing of the AtmoSwing (http://www.atmoswing.org/) results in R

You can use devtools to install the package:

```r
install.packages('devtools')
library('devtools')
install_github(repo='atmoswing/R-toolbox',dependencies=T)
library(atmoswing)
?atmoswing
```

## Development

To generate the docs, do:

```r
roxygen2::roxygenise()
```

To start the tests, do:

```r
devtools::test()
```
