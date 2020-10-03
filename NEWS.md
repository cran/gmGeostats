# gmGeostats 0.10.7

* (2020-09-27) updated references in generalised diagonalisations; improved/reimplemented `as.logratioVariogram()` and `as.logratioVariogramAnisotropy()` methods; minor edits at `plot.logratioVariogramAnisotropy()` and `variogramModelPlot.logratioVariogram()`
* (2020-09-24) openMP pragmas enclosed, to improve portability 
* (2020-09-18) conversion of variogram models of class CompLinModCoReg -> LMCAnisCompo -> variogramModel(List)
* (2020-09-18) `gsi.gstatCokriging2compo()` and `gsi.gstatCokriging2rmult()` now work with systems of 3 components resp. 2 variables (systems of 2 components or 1 variable still fail)

# gmGeostats 0.10.6

* (2020-09-17) we are in CRAN! https://cran.r-project.org/web/packages/gmGeostats/