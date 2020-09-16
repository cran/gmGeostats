

#### spatial method specifications --------------
# S3 -> S4 classes
# cat("creating spatial method parameter classes\n")

#' @include gmAnisotropy.R
#' @include gmValidationStrategy.R 
#' @include variograms.R
setOldClass("gmKrigingNeighbourhood")
setOldClass("gmDirectSamplingParameters")
setOldClass("gmTurningBands")
setOldClass("gmSequentialSimulation")
setOldClass("gmCholeskyDecomposition")
setOldClass("NfoldCrossValidation")
setOldClass("LeaveOneOut")

# abstract classes
# cat("creating spatial method parameter classes: superclass creation\n")

#' @title Neighbourhood description
#' @description abstract class, containing any specification of a spatial neighbourhood
#' @export
setClassUnion(name="gmNeighbourhoodSpecification", members=c("gmKrigingNeighbourhood","NULL"))

#' @title Validation strategy description
#' @description abstract class, containing any specification of a validation strategy for spatial models
#' @export
setClassUnion(name="gmValidationStrategy", 
              members=c("NULL",
                        "LeaveOneOut", 
                        "NfoldCrossValidation"))

#' @title parameters for Multiple-Point Statistics methods
#' @description abstract class, containing any parameter specification of a spatial multipoint algorithm 
#' @export
setClassUnion(name="gmMPSParameters", 
              members=c("gmDirectSamplingParameters","NULL"))


#' @title parameters for Gaussian Simulation methods
#' @description abstract class, containing any parameter specification of a spatial simulation algorithm
#' exploiting a Gaussian two-point model structure 
#' @export
setClassUnion(name="gmGaussianSimulationAlgorithm", 
              members=c("gmSequentialSimulation", 
                        "gmTurningBands",
                        "gmCholeskyDecomposition",
                        "NULL") )

#' @title Parameter specification for a spatial simulation algorithm
#' @description abstract class, containing any parameter specification for a spatial simulation algorithm
#' @export
setClassUnion(name="gmSimulationAlgorithm", 
              members=c("gmGaussianSimulationAlgorithm", 
                        "gmMPSParameters"))

#' @title parameters for Spatial Gaussian methods of any kind
#' @description abstract class, containing any parameter specification for a spatial algorithm 
#' for interpolation, simulation or validation making use of Gaussian assumptions
#' @export
setClassUnion(name="gmGaussianMethodParameters", 
              members=c("gmSequentialSimulation", 
                        "gmKrigingNeighbourhood",
                        "gmValidationStrategy"))


#' @title Parameter specification for any spatial method
#' @description abstract class, containing any parameter specification for 
#' any spatial method. Members of this class are [gmNeighbourhoodSpecification-class]
#' [gmMPSParameters-class] and [gmValidationStrategy-class].
#' 
#' @export
setClassUnion(name="gmSpatialMethodParameters", 
              members=c("NULL",
                        "gmNeighbourhoodSpecification",
                        "gmMPSParameters",
                        "gmValidationStrategy")
)





#' @title MPS training image class
#' @description abstract class, containing any specification of a multiple-point
#' training image
#' @export
#' @importClassesFrom sp SpatialGridDataFrame SpatialGrid GridTopology
#' @importClassesFrom sp SpatialPixelsDataFrame SpatialPixels
#' @importClassesFrom sp SpatialPointsDataFrame SpatialPoints
setClassUnion(name="gmTrainingImage", 
              members=c("SpatialGridDataFrame", 
                        "SpatialPixelsDataFrame")
)


#' @title General description of a spatial model
#' @description abstract class, containing any specification of an unconditional
#' spatial model 
#' @export
setClassUnion(name="gmUnconditionalSpatialModel", 
              members=c("NULL",
                        "gmGaussianModel",
                        "gmTrainingImage")
)

