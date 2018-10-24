
# nipals 0.4 - Oct 2018

Minor fix: add row/column names to fitted matrix.

New: `startcol` argument can now be a function.

Slight change to automatic start column selection. By default, the start column is the column with the largest sum of absolute values. (Formerly was largest variance.)

# nipals 0.3 - Nov 2017

The `nipals` function was split off from the `gge` package, extensively optimized and compared with implementations in other packages.

## Todo

Find published example of NIPALS with missing data. (Only found 1)
