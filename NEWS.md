
# nipals 0.6 - Jan 2010

* New (experimental) function `empca()` to calculate principal components via weighted EM PCA.

* Calculating fitted values with `ncomp=1` now works (#2)

* `nipals()` now returns `center` and `scale` in the fitted model object, which makes it easier to get predictions from the model (#1).

# nipals 0.5 - Oct 2018

* Row/column names added to fitted matrix.

* New `startcol` argument can now be a function which is applied to every column.

* There is a slight change to automatic start column selection. By default, the start column is the column with the largest sum of absolute values. (Formerly, the start column was the column with the largest variance, but this does not make sense when columns are scaled.)

# nipals 0.4 - Feb 2018

# nipals 0.3 - Nov 2017

* The `nipals()` function was split off from the `gge` package, extensively optimized and compared with implementations in other packages.
