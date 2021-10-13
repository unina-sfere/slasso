## Test environments
* local R installation, R 4.0.5
* GitHub Actions (ubuntu):   devel, release, oldrel-1
* GitHub Actions (windows): release
* Github Actions (macOS): release

## R CMD check results

0 errors | 0 warnings | 2 note

* This is a new release.
* Check: package dependencies, Result: NOTE Package in Depends/Imports which should probably only be in LinkingTo: 'RcppArmadillo'. RccpArmadillo in the Imports field is needed to succesfully call cxxfunctionplus2 in the zzz.R file for multiple architecture installation.
