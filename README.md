# nipals <img src="figure/nipals_logo_150.png" align="right" />

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/nipals)](https://cran.r-project.org/package=nipals)
[![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/nipals)](https://cranlogs.r-pkg.org/badges/nipals)
[![Research software impact](http://depsy.org/api/package/cran/nipals/badge.svg)](http://depsy.org/package/r/nipals)

Key features:
  
* Allows missing values in the data.
* Uses Gram-Schmidt to ensure orthogonal principal components.
* Carefully optimized for speed.
* Flexible options.
* Unit tests.

## Installation

```R
# Install the released version from CRAN:
install.packages("nipals")

# Install the development version from GitHub:
install.packages("devtools")
devtools::install_github("kwstat/nipals")
```

## Usage

Vignettes:
  [Introduction to the nipals package](https://rawgit.com/kwstat/nipals/master/vignettes/nipals_examples.html)

```R
require(nipals)
data(uscrime, package = "nipals")

```
