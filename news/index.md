# Changelog

## nipals 1.2 (2026-09-08)

CRAN release: 2026-09-09

- Fix CRAN check notes.

- Use Air formatter

- Reviewed by Claude Opus.

## nipals 1.0 (2024-12-02)

CRAN release: 2024-12-02

- Add
  [`avg_angular_distance()`](https://kwstat.github.io/nipals/reference/avg_angular_distance.md)
  function.

- Add more tests.

- Switch to MIT license.

- Documentation pages now created via Github Actions.

- Clarify that
  [`nipals()`](https://kwstat.github.io/nipals/reference/nipals.md)
  returns `T*Lambda*P'`.

## nipals 0.8 (2021-08-01)

CRAN release: 2021-09-15

- Fix issue ([\#5](https://github.com/kwstat/nipals/issues/5)). A
  floating-point error created a dot product outside \[-1,1\].

## nipals 0.7 (2020-01-24)

CRAN release: 2020-01-24

- Fix ‘additional issues’ reported by CRAN
  ([\#3](https://github.com/kwstat/nipals/issues/3)).

## nipals 0.6 (2020-01-01)

- Calculating xhat with `ncomp=1` now works
  ([\#2](https://github.com/kwstat/nipals/issues/2)).

- [`nipals()`](https://kwstat.github.io/nipals/reference/nipals.md) now
  returns `center` and `scale` in the fitted model object, which makes
  it easier to get predictions from the model
  ([\#1](https://github.com/kwstat/nipals/issues/1)).

- New (experimental) function
  [`empca()`](https://kwstat.github.io/nipals/reference/empca.md) to
  calculate principal components via weighted EM PCA.

- Row/column names added to fitted matrix.

## nipals 0.5 (2018-10-24)

CRAN release: 2018-10-24

- Row/column names added to fitted matrix.

- New `startcol` argument can now be a function which is applied to
  every column.

- There is a slight change to automatic start column selection. By
  default, the start column is the column with the largest sum of
  absolute values. (Formerly, the start column was the column with the
  largest variance, but this does not make sense when columns are
  scaled.)

## nipals 0.4 (2018-02-11)

CRAN release: 2018-02-11

## nipals 0.3 (2017-11-14)

CRAN release: 2017-11-14

- The [`nipals()`](https://kwstat.github.io/nipals/reference/nipals.md)
  function was split off from the `gge` package, extensively optimized
  and compared with implementations in other packages.

- First release to CRAN.
