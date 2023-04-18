# gmGeostats 0.11.3

* (2023-04-13) batch of hidden methods exported to the namespace (mostly for `as.*()` functions converting between different representations of empirical variograms and of geostatistical models)
* (2023-04-04) C routines called as symbols and not name strings; non-ASCII characters removed

# gmGeostats 0.11.2

* (2022-11-22) flag '-lstdc++' removed from makevars.
* (2022-11-15) minor bugs in accuracy calculations corrected. References improved.
* (2022-08-15) bug in turning bands for spherical variogram corrected.


# gmGeostats 0.11.1

* (2022-05-03) minor changes to adapt to more stringent standards for class determination, for variable definition in C code, and for calls to FORTRAN routines from C code involving string arguments.

# gmGeostats 0.11.0-9002

* (2021-12-14) bugs in turning bands corrected: getUnitVec was producing a wrong sequence of directions, and bands for exponential and Gaussian variograms did not correct for the difference between parametric and effective range in the right way; 
* (2021-12-08) anisotropy objects complete, restructured and everywhere correctly used: if `class(A)=="AnisotropyScaling"` and `class(M)=="AnisotropyRangeMatrix"` (the two possible classes), then `A %*% t(A) == solve(M)`; with them, `u = sqrt(h %*% solve(M, h))` can be used in a scalar variogram (i.e. `M` is a matrix of ranges and orientations), and `v = t(A) %*% h` an isotropic variogram (i.e. `u=sqrt(sum(v^2))`,  and `A` is a scaling matrix that makes the space isotropic)
* (2021-12-08) bugs corrected in `gsi.EVario3D()`
 
# gmGeostats 0.11.0-9001

* (2021-11-05) bugs corrected in internal 2D empirical variogram function `gsi.EVario2D()` . First version of 3D internal variogram function `gsi.EVario3D()` available. Usage of these functions is only for specialists foreseen. In the near future a user-friendlier wrapper will be provided.

# gmGeostats 0.11.0-9000

* (2021-10-20) section on the vignette "register_new_layer_datatype.Rmd" about the definition and registration of empirical covariance for circular data.

# gmGeostats 0.11.0

* (2021-10-17) dependence from randomFields eliminated 
* (2021-10-17) consolidated management of `predict()` methods for S3 and S4 objects, with the help of the internal function `Predict()` (not recommended to use)

# gmGeostats 0.10.9-9001

* (2021-09-30) minor bug corrected on the way fit_lmc passes its arguments to surrogate gstat and gmGeostats functions
* (2021-10-01) management of `predict()` methods for S3 and S4 objects, with the help of the internal function `Predict()` (not recommended to use)

# gmGeostats 0.10.9

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
