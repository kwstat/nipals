# nipals <img src="man/figures/logo.png" align="right" />

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/nipals)](https://cran.r-project.org/package=nipals)
[![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/nipals)](https://cranlogs.r-pkg.org/badges/nipals)

Homepage: https://kwstat.github.io/nipals

Repository: https://github.com/kwstat/nipals

The `nipals` package provides two functions to perform Principal Components Analysis of a matrix:
1. The `nipals` function uses Non-linear Iterative Partial Least Squares.
2. The `empca` function uses EM PCA.

NIPALS has been implemented several times in R packages. EMPCA has previously appeared in python, but is available in R here for the first time.
This package strives to have the best (fast and accurate) R implementations.

The `empca()` function should be considered **experimental**.  There is a problem using `empca()` on matrices that are both (1) non-full rank (2) have missing values.

## Key features

* Missing values are allowed.
* Uses Gram-Schmidt to ensure orthogonal principal components.
* Carefully optimized for speed (`nipals` only, not `empca`).
* Flexible options.
* Vignettes and unit tests.
* Weights are allowed (`empca` only).

## Installation

```R
# Install the released version from CRAN:
install.packages("nipals")

# Install the development version from GitHub:
install.packages("devtools")
devtools::install_github("kwstat/nipals")
```

## Usage

```R
require(nipals)
data(uscrime, package = "nipals")
dat <- uscrime
dat <- as.matrix(dat[ , -1])

# Gram-Schmidt corrected NIPALS
m3 <- nipals(dat)

# Show that the Principal Components are orthogonal
round(crossprod(m3$loadings),3)
##      PC1 PC2 PC3 PC4 PC5 PC6 PC7
##  PC1   1   0   0   0   0   0   0
##  PC2   0   1   0   0   0   0   0
##  PC3   0   0   1   0   0   0   0
##  PC4   0   0   0   1   0   0   0
##  PC5   0   0   0   0   1   0   0
##  PC6   0   0   0   0   0   1   0
##  PC7   0   0   0   0   0   0   1

round(m3$loadings,3)
##              PC1    PC2    PC3    PC4    PC5    PC6    PC7
##  murder    0.296 -0.623  0.178 -0.241  0.540 -0.265 -0.270
##  rape      0.432 -0.171 -0.244  0.060  0.200  0.769  0.299
##  robbery   0.397  0.044  0.496 -0.558 -0.519  0.120  0.005
##  assault   0.399 -0.353 -0.063  0.629 -0.502 -0.166 -0.191
##  burglary  0.440  0.204 -0.211 -0.057  0.095 -0.540  0.645
##  larceny   0.358  0.401 -0.541 -0.231  0.023 -0.037 -0.602
##  autotheft 0.295  0.504  0.567  0.418  0.372  0.054 -0.148

m4 <- empca(dat)

# Show that the Principal Components are orthogonal
round(crossprod(m4$loadings),3)
##     PC1 PC2 PC3 PC4 PC5 PC6 PC7
## PC1   1   0   0   0   0   0   0
## PC2   0   1   0   0   0   0   0
## PC3   0   0   1   0   0   0   0
## PC4   0   0   0   1   0   0   0
## PC5   0   0   0   0   1   0   0
## PC6   0   0   0   0   0   1   0
## PC7   0   0   0   0   0   0   1

round(m4$loadings,3)
##             PC1    PC2    PC3    PC4    PC5    PC6    PC7
## murder    0.301 -0.642  0.168 -0.433  0.447  0.062  0.278
## rape      0.431 -0.146 -0.266  0.106 -0.018 -0.799 -0.269
## robbery   0.399  0.023  0.488 -0.306 -0.712  0.028 -0.012
## assault   0.396 -0.336 -0.099  0.734 -0.161  0.349  0.181
## burglary  0.441  0.197 -0.191 -0.179  0.229  0.459 -0.659
## larceny   0.353  0.430 -0.515 -0.245 -0.044  0.050  0.601
## autotheft 0.295  0.478  0.593  0.275  0.461 -0.150  0.150

```

## See also

A python version of this package can be found at <https://pypi.org/project/nipals/>.
