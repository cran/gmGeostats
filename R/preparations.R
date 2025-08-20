### generic S4 functions --------------
#' @export
#' @rdname predict_gmSpatialModel
setGeneric("Predict", function(object, newdata, pars, ...){
  standardGeneric("Predict")
})

#' @export
#' @rdname predict_gmSpatialModel
setGeneric("predict", function(object,...) standardGeneric("predict"))  



#' Empirical variogram or covariance function 
#'
#' Generic function to compute the empirical variogram or covariance function from actual data
#'
#' @param object spatial data container, with special methods according to the class  
#' @param ... further parameters for generic functionality
#'
#' @return An empirical variogram for the provided data. 
#' @importFrom gstat variogram
#' @export
#' @family gmEVario functions
setGeneric("variogram", function(object,...) standardGeneric("variogram"))  


# @include gmAnisotropy.R
# @include gmValidationStrategy.R 
# @include variograms.R
setOldClass("gmKrigingNeighbourhood")
setOldClass("gmDirectSamplingParameters")
setOldClass("gmTurningBands")
setOldClass("gmSequentialSimulation")
setOldClass("gmCholeskyDecomposition")
setOldClass("LeaveOneOut")
setOldClass("NfoldCrossValidation")



# @include gstatCompatibility.R
# @include compositionsCompatibility.R
setOldClass("gmCgram")
setOldClass("LMCAnisCompo")
setOldClass("variogramModelList")
setOldClass("variogramModel")
setOldClass("genDiag")


# @include gstatCompatibility.R
# @include compositionsCompatibility.R
setOldClass("gmEVario")



#' Logratio variogram of a compositional data
#' 
#' gmGeostats reimplementation of the compositions::logratioVariogram function
#'
#' @inheritParams logratioVariogram.default 
#'
#' @return a "logratioVariogram" object, or a "logratioVariogramAnisotropy" object
#' if you provide more than one `azimuth`. See [logratioVariogram()] for details and 
#' examples.
#' @export
setOldClass("logratioVariogram")
setOldClass("logratioVariogramAnisotropy")
setOldClass("gstatVariogram")


