# nipals 0.9 - unpublished

* Add `avg_angular_distance()` function.
* Add more tests.
* Switch to MIT license.

# nipals 0.8 - Sep 2021

* Fix issue (#5). A floating-point error created a dot product outside [-1,1].

# nipals 0.7 - Jan 2020

* Fix 'additional issues' reported by CRAN (#3).

# nipals 0.6 - Jan 2020

* Calculating xhat with `ncomp=1` now works (#2).

* `nipals()` now returns `center` and `scale` in the fitted model object, which makes it easier to get predictions from the model (#1).

* New (experimental) function `empca()` to calculate principal components via weighted EM PCA.

* Row/column names added to fitted matrix.

# nipals 0.5 - Oct 2018

* Row/column names added to fitted matrix.

* New `startcol` argument can now be a function which is applied to every column.

* There is a slight change to automatic start column selection. By default, the start column is the column with the largest sum of absolute values. (Formerly, the start column was the column with the largest variance, but this does not make sense when columns are scaled.)

# nipals 0.4 - Feb 2018


# nipals 0.3 - Nov 2017

* The `nipals()` function was split off from the `gge` package, extensively optimized and compared with implementations in other packages.
* First release to CRAN.
