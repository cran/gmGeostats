#' @import methods compositions RColorBrewer
NULL

# 
# #### variogram functions -------------
#  
# ## theoretical structural functions
## S3 -> S4 classes  
# # cat("creating variogram model classes\n")
# setOldClass("gmCgram")
# setOldClass("LMCAnisCompo")
# setOldClass("variogramModelList")
# setOldClass("variogramModel")
# 
# # abstract classes
# #' @title Structural function model specification
# #' @description Abstract class, containing any specification of a variogram (or covariance) model
# #' @export
# setClassUnion(name="ModelStructuralFunctionSpecification", 
#               members=c("NULL","gmCgram", "LMCAnisCompo", "variogramModelList", "variogramModel"))
# 
# 
# 
# 
# 
# #### container class --------------
# #' An S4 class to represent a Gaussian random field specification
# #'
# #' @slot structure ModelStructuralFunctionSpecification. Variogram or 
# #' (generalised) covariance function specification, typically an object
# #' obtained from a call to functions such as \code{\link{setCgram}},
# #' \code{\link{LMCAnisCompo}} or \code{gstat::vgm}.
# #' @slot formula formula specifying the structure
# #' of dependence of the mean of the random field w.r.to spatial coordinates 
# #' and/or covariables; typically it will have no left-hand-side term; 
# #' @slot beta numeric, a vector with as many coefficients as terms the formula
# #' above requires for a full specification of the trend; if unknown, these can
# #' be NAs, as many as needed. 
# #'
# #' @return A object with the slots populated as given
# #' @export
# #' @seealso [gmSpatialModel-class], and the `make.gm*` functions referenced there
# setClass("gmGaussianModel", 
#          slots = list(structure = "ModelStructuralFunctionSpecification",
#                       formula="formula",
#                       beta = "structure")
# )
# 
# 
# ## empirical structural functions
# # S3 -> S4 classes
# # cat("creating empirical variogram classes\n")
# setOldClass("gmEVario")
# setOldClass("logratioVariogram")
# setOldClass("logratioVariogramAnisotropy")
# setOldClass("gstatVariogram")
# 
# 
# # abstract classes
# #' @title Empirical structural function specification
# #' @description Abstract class, containing any specification of an empirical variogram (or covariance function, or variations)
# #' @export
# setClassUnion(name="EmpiricalStructuralFunctionSpecification", members=c("NULL","gmEVario", "logratioVariogram", "logratioVariogramAnisotropy", "gstatVariogram"))
# 
# 
# 
# #' #### spatial method specifications --------------
# #' # S3 -> S4 classes
# #' # cat("creating spatial method parameter classes\n")
# #' 
# #' setOldClass("gmKrigingNeighbourhood")
# #' setOldClass("gmDirectSamplingParameters")
# #' setOldClass("gmTurningBands")
# #' setOldClass("gmSequentialSimulation")
# #' setOldClass("gmCholeskyDecomposition")
# #' setOldClass("NfoldCrossValidation")
# #' setOldClass("LeaveOneOut")
# #' 
# #' # abstract classes
# #' # cat("creating spatial method parameter classes: superclass creation\n")
# #' 
# # #' @title Neighbourhood description
# # #' @description abstract class, containing any specification of a spatial neighbourhood
# # #' @export
# #' setClassUnion(name="gmNeighbourhoodSpecification", members=c("gmKrigingNeighbourhood","NULL"))
# #' 
# # #' @title Validation strategy description
# # #' @description abstract class, containing any specification of a validation strategy for spatial models
# # #' @export
# #' setClassUnion(name="gmValidationStrategy", 
# #'               members=c("NULL",
# #'                         "LeaveOneOut", 
# #'                         "NfoldCrossValidation"))
# #' 
# # #' @title parameters for Multiple-Point Statistics methods
# # #' @description abstract class, containing any parameter specification of a spatial multipoint algorithm 
# # #' @export
# #' setClassUnion(name="gmMPSParameters", 
# #'               members=c("gmDirectSamplingParameters","NULL"))
# #' 
# #' 
# # #' @title parameters for Gaussian Simulation methods
# # #' @description abstract class, containing any parameter specification of a spatial simulation algorithm
# # #' exploiting a Gaussian two-point model structure 
# # #' @export
# #' setClassUnion(name="gmGaussianSimulationAlgorithm", 
# #'               members=c("gmSequentialSimulation", 
# #'                         "gmTurningBands",
# #'                         "gmCholeskyDecomposition",
# #'                         "NULL") )
# #' 
# # #' @title Parameter specification for a spatial simulation algorithm
# # #' @description abstract class, containing any parameter specification for a spatial simulation algorithm
# # #' @export
# #' setClassUnion(name="gmSimulationAlgorithm", 
# #'               members=c("gmGaussianSimulationAlgorithm", 
# #'                         "gmMPSParameters"))
# #' 
## #' @title parameters for Spatial Gaussian methods of any kind
## #' @description abstract class, containing any parameter specification for a spatial algorithm 
## #' for interpolation, simulation or validation making use of Gaussian assumptions
## #' @export
## setClassUnion(name="gmGaussianMethodParameters", 
##               members=c("gmSequentialSimulation", 
##                         "gmKrigingNeighbourhood",
##                         "gmValidationStrategy"))
##
## 
## #' @title Parameter specification for any spatial method
## #' @description abstract class, containing any parameter specification for any spatial method
## #' @export
## setClassUnion(name="gmSpatialMethodParameters", 
##               members=c("NULL",
##                         "gmNeighbourhoodSpecification",
##                         "gmMPSParameters",
##                         "gmValidationStrategy")
## )
# 
# 
# #' @title MPS training image class
# #' @description abstract class, containing any specification of a multiple-point
# #' training image
# #' @export
#' setClassUnion(name="gmTrainingImage", 
#'               members=c("SpatialGridDataFrame", 
#'                         "SpatialPixelsDataFrame")
#' )
#' 
#' 
# # @title General description of a spatial model
# # @description abstract class, containing any specification of an unconditional
# #' spatial model 
# #' @export
# setClassUnion(name="gmUnconditionalSpatialModel", 
#               members=c("NULL",
#                         "gmGaussianModel",
#                         "gmTrainingImage")
# )
# 

#' @title General description of a spatial data container
#' @description abstract class, containing any specification of a spatial data container
#' @export
setClassUnion(name="gmSpatialDataContainer", 
              members=c("NULL",
                        "SpatialPointsDataFrame",
                        "SpatialPixelsDataFrame",
                        "SpatialGridDataFrame")
)




### data containers -----------------
# S3 -> S4 classes
#cat("creating complex data container classes\n")
#setOldClass("gmMultiDataFrame", S4Class="data.frame")
setOldClass(c("DataFrameStack", S4Class="data.frame"))  

# abstract classes
#' @title Superclass for grid or nothing
#' @description Superclass for slots containing a grid topology or being empty
#' @export
setClassUnion(name="GridOrNothing", members = c("NULL", "GridTopology"))



.onAttach <- function(libname, pkgname, ...){
  utils::data("variogramModels")
  
  ## package startup message  
  packageStartupMessage("Welcome to 'gmGeostats', a package for multivariate geostatistical analysis.\n Note: use 'fit_lmc' instead of fit.lmc")
}


# vg.Gau <- vg.Gauss <- vg.gauss <- 0
# vg.Sph <- vg.Spherical <- vg.sph <- 1
# vg.Exp <- vg.Exponential <- vg.exp <- 2
# gsi.validModels <- 0:2


.onLoad <- function(libname, pkgname){

  ## set package options ---- 
    # grid organisation
    gridOrder = list(refpoint="topleft", cycle=1:2)
    scaleClasses = list(continuous=c("acomp","aplus", "rcomp", "rplus", "rmult"), discrete=c("ccomp", "factor"))
    o = list(gridOrder=gridOrder, scaleClasses=scaleClasses)
    options(gmGeostats=o)
  
  gsi.validModels <- 0:2

  ## set up generic functionality ---
  # if(!exists("fit.lmc") | !isGeneric("fit.lmc")) fit.lmc <- function(v, ...) UseMethod("fit.lmc", v)
  
  invisible()  
}

.onUnload <- function(libpath){
  library.dynam.unload("gmGeostats", libpath)
  ## remove package options ---- 
  options(gmGeostats=NULL)
  
  invisible()  
}
