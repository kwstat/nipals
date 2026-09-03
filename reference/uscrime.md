# U.S. Crime rates per 100,00 people

U.S. Crime rates per 100,00 people for 7 categories in each of the 50
U.S. states in 1977.

## Usage

``` r
uscrime
```

## Format

A data frame with 50 observations on the following 8 variables.

- state:

  U.S. state

- murder:

  murders

- rape:

  rapes

- robbery:

  robbery

- assault:

  assault

- burglary:

  burglary

- larceny:

  larceny

- autotheft:

  automobile thefts

## Source

Documentation Example 3 for PROC HPPRINCOMP.
http://documentation.sas.com/api/docsets/stathpug/14.2/content/stathpug_code_hppriex3.htm?locale=en

## Details

There are two missing values.

## References

SAS/STAT User's Guide: High-Performance Procedures. The HPPRINCOMP
Procedure.
http://support.sas.com/documentation/cdl/en/stathpug/67524/HTML/default/viewer.htm#stathpug_hpprincomp_toc.htm

## Examples

``` r

library(nipals)
head(uscrime)
#>        state murder rape robbery assault burglary larceny autotheft
#> 1    Alabama   14.2 25.2    96.8   278.3   1135.5  1881.9     280.7
#> 2     Alaska   10.8 51.6    96.8   284.0   1331.7  3369.8     753.3
#> 3    Arizona    9.5 34.2   138.2   312.3   2346.1  4467.4     439.5
#> 4   Arkansas    8.8 27.6    83.2   203.4    972.6  1862.1     183.4
#> 5 California   11.5 49.4   287.0   358.0   2139.4  3499.8     663.5
#> 6   Colorado    6.3 42.0   170.7   292.9   1935.2  3903.2     477.1

# SAS deletes rows with missing values
dat <- uscrime[complete.cases(uscrime), ]
dat <- as.matrix(dat[ , -1])
m1 <- nipals(dat) # complete-data method

# Traditional NIPALS with missing data  
dat <- uscrime
dat <- as.matrix(dat[ , -1])
m2 <- nipals(dat, gramschmidt=FALSE) # missing 
round(crossprod(m2$loadings),3) # Prin Comps not quite orthogonal
#>        PC1    PC2    PC3    PC4    PC5    PC6    PC7
#> PC1  1.000  0.002  0.002 -0.001 -0.001 -0.002  0.005
#> PC2  0.002  1.000  0.001 -0.001  0.000 -0.006  0.004
#> PC3  0.002  0.001  1.000  0.000  0.000 -0.001  0.000
#> PC4 -0.001 -0.001  0.000  1.000  0.000  0.001  0.002
#> PC5 -0.001  0.000  0.000  0.000  1.000 -0.002 -0.002
#> PC6 -0.002 -0.006 -0.001  0.001 -0.002  1.000  0.005
#> PC7  0.005  0.004  0.000  0.002 -0.002  0.005  1.000
  
# Gram-Schmidt corrected NIPALS
m3 <- nipals(dat, gramschmidt=TRUE) # TRUE is default
round(crossprod(m3$loadings),3) # Prin Comps are orthogonal
#>     PC1 PC2 PC3 PC4 PC5 PC6 PC7
#> PC1   1   0   0   0   0   0   0
#> PC2   0   1   0   0   0   0   0
#> PC3   0   0   1   0   0   0   0
#> PC4   0   0   0   1   0   0   0
#> PC5   0   0   0   0   1   0   0
#> PC6   0   0   0   0   0   1   0
#> PC7   0   0   0   0   0   0   1
```
