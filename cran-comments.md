## Test environments
* Mac OS 15.3.1, R 4.4.1
* Windows 10 Pro, R 4.5.1
* win-builder (devel and release)

## R CMD check results

** Mac OS 15.3.1, R 4.5.0
0 errors | 0 warnings | 1 note | 1 info

* checking CRAN incoming feasibility ... [3s/10s] NOTE
Maintainer: ‘Silvia Castiglione <silvia.castiglione@unina.it>’
Size of tarball: 6548969 bytes

* checking installed package size ... INFO
  installed size is  9.0Mb
  sub-directories of 1Mb or more:
    doc   7.1Mb
    
-This is just memory allocation to vignette files.


* checking for future file timestamps ... NOTE
unable to verify current time

-This is an erratic note impossible to resolve locally.



** Windows 10 Pro, R 4.5.1
0 errors | 0 warnings | 1 note

* checking for future file timestamps ... NOTE
unable to verify current time

-This is an erratic note impossible to resolve locally.



** win-builder (devel and release)
** release
0 errors | 0 warnings | 0 note
STATUS OK


** devel
0 errors | 0 warnings | 0 note
STATUS OK

## Reverse dependencies

We checked 2 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

