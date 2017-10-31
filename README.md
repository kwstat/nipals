# nipals <img src="figure/nipals_logo_150.png" align="right" />

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/nipals)](https://cran.r-project.org/package=nipals)
[![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/nipals)](https://cranlogs.r-pkg.org/badges/nipals)
[![Research software impact](http://depsy.org/api/package/cran/nipals/badge.svg)](http://depsy.org/package/r/nipals)

The 'nipals' package provides a single function to perform Principal Components Analysis of a matrix using Non-linear Iterative Partial Least Squares. NIPALS has been implemented several times in R packages.  This package strives to be the best implementation.

Key features:
  
* Missing values are allowed.
* Uses Gram-Schmidt to ensure orthogonal principal components.
* Carefully optimized for speed.
* Flexible options.
* Vignettes and unit tests.

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

  [NIPALS algorithm](https://rawgit.com/kwstat/nipals/master/vignettes/nipals_algorithm.html)

  [Comparing NIPALS functions in R](https://rawgit.com/kwstat/nipals/master/vignettes/nipals_comparisons.html)

  [NIPALS optimization notes](https://rawgit.com/kwstat/nipals/master/vignettes/nipals_optimization.html)

```R
require(nipals)
data(uscrime, package = "nipals")

```
