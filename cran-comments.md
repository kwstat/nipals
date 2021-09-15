# nipals 0.8

This fixes an error identified by CRAN on Mac M1.  See: https://github.com/kwstat/nipals/issues/5

* Rhub: Mac M1

# nipals 0.7

This is a resubmission. Thanks to Uwe Ligges for patiently explaining the 'Additional issues' caused by --noLD.  I fixed the problem and verified the fix via the Rhub platform "debian-gcc-devel-nold".

## Test environments

* local: Windows 10, R 3.6.2
* Rhub: Debian Linux, R-devel, GCC, no long double
* Win-builder: R-release, version 3.6.2 (2019-12-12)
* Win-builder: R-devel ATC

## Rcmd check results

No errors or warnings.  One note:

  Possibly mis-spelled words in DESCRIPTION:
    EMPCA (2:63)

This abbreviation is correctly spelled.

## Downstream dependencies

gge: Checked, no problems.


# nipals 0.4

## Test environments

* local: Windows 7, R 3.5.1
* win-builder: R-release 3.5.1
* win-builder: R-devel

## Rcmd check results

No errors, warnings, or notes.
  
## Downstream dependencies

gge: No problems.

Downloading 1 source packages: gge
Checking packages --------------------------------------------------------------------------
Checking 1 packages: gge
Package library: C:/Temp/Rtmp2nviSi/revdepfb05b6390a, C:/Temp/Rtmp2nviSi/R-lib, C:/Program Files/R/R-3.5.1/library
Checked gge: 0 errors | 0 warnings | 0 notes

# nipals 0.2

## Reviewer comments

Thanks, please add a reference for the method in the 'Description' field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
with no space after 'doi:', 'arXiv:' and angle brackets for auto-linking.
