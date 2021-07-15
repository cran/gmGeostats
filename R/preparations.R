### generic S4 functions --------------
#' @export
#' @rdname predict_gmSpatialModel
setGeneric("Predict", function(object, newdata, pars, ...){
  standardGeneric("Predict")
})

#' @export
#' @rdname predict_gmSpatialModel
setGeneric("predict", function(object,...) standardGeneric("predict"))  



#' @export
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
setOldClass("logratioVariogram")
setOldClass("logratioVariogramAnisotropy")
setOldClass("gstatVariogram")


