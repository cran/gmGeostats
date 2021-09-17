# gmGeostats 0.10-9

(release of polished 0.10.8.900x dev versions)
* minor correction on C source code for compatibility with clang13

# gmGeostats 0.10.8-9001

* (2021-07-22) conversion methods between variogram models: `as.gmCgram()` methods for "variogramModel" and "variogramModelList" objects (from package "gstat")

# gmGeostats 0.10.8-9000

* (2021-07-15) documented abstract union classes, specifying required methods for the member classes

# gmGeostats 0.10.8

(release of polished 0.10.7.900x dev versions)

# gmGeostats 0.10.7.9006

* new vignette "How to register new layer datatypes", with a howto adapt gmGeostats to deal with your very special data type (with an example for regionalized circular data)
* (2021-07-06) new option to extract univariate accuracy() from multivariate kriging results
* NEWS.md and vignettes/gmGeostats.Rmd polished
* "predict"-methods for objects of class gmSpatialModel now programmed with multiple dispatching, to enable further extensions

# gmGeostats 0.10.7.9005

* (2021-06-30) bugs in xvErrorMeasures() for simulated data corrected
* (2021-06-30) bugs in accuracy() and precision() corrected

# gmGeostats 0.10.7.9004

* (2021-06-29) bugs in accuracy() and precision() corrected
* (2021-06-23) bug in validate.NfoldCrossValidation() corrected

# gmGeostats 0.10.7.9003

* (2021-04-08) exported turning bands for direct usage
* (2020-10-01) bug in turning bands spherical variogram corrected

#  gmGeostats 0.10.7

* (2020-09-27) updated references in generalised diagonalisations; improved/reimplemented `as.logratioVariogram()` and `as.logratioVariogramAnisotropy()` methods; minor edits at `plot.logratioVariogramAnisotropy()` and `variogramModelPlot.logratioVariogram()`
* (2020-09-24) openMP pragmas enclosed, to improve portability 
* (2020-09-18) conversion of variogram models of class CompLinModCoReg -> LMCAnisCompo -> variogramModel(List)
* (2020-09-18) `gsi.gstatCokriging2compo()` and `gsi.gstatCokriging2rmult()` now work with systems of 3 components resp. 2 variables (systems of 2 components or 1 variable still fail)

# gmGeostats 0.10.6

* (2020-09-17) we are in CRAN! https://cran.r-project.org/web/packages/gmGeostats/
