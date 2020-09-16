##### KrigingNeighbourhood --------
#' Create a parameter set of local for neighbourhood specification.
#'
#' Create a parameter set describing a kriging neighbourhood (local or global) for 
#' cokriging and cokriging based simulation. This heavily relies on the definitions of 
#' [gstat::gstat()]. All parameters are optional, as their default amounts to a global 
#' neihghbourhood.
#'
#' @param nmax maximum number of data points per cokriging system
#' @param nmin minimum number of data points per cokriging system
#' @param omax maximum number of data points per cokriging system per quadrant/octant
#' @param maxdist maximum radius of the search neighborhood
#' @param force logical; if less than `nmin` points are found inside `maxdist` radius, 
#' keep searching. 
#' @param anisotropy currently ignored; in the future, argument to specify anisotropic search areas.
#' @param ... further arguments, currently ignored
#'
#' @return an S3-list of class "gmKrigingNeighbourhood" containing the six elements given as arguments 
#' to the function. This is just a compact way to provide further functions such as [predict.gmSpatialModel()]
#' with appropriate triggers for choosing a prediction method or another, in this case for triggering
#' cokriging (if alone) or eventually sequential simulation (see [SequentialSimulation()]).
#' @export
#'
#' @examples
#' data("jura", package="gstat")
#' X = jura.pred[,1:2]
#' summary(X)
#' Zc = jura.pred[,7:10]
#' ng_global = KrigingNeighbourhood()
#' ng_local = KrigingNeighbourhood(maxdist=1, nmin=4, 
#'                                 omax=5, force=TRUE)
#' ng_local
#' ng_global
#' make.gmCompositionalGaussianSpatialModel(data = Zc, coords = X,
#'                                          V = "alr", ng = ng_local)
KrigingNeighbourhood <- function(nmax=Inf, nmin=0, omax=0, maxdist=Inf, force=FALSE, anisotropy=NULL, ...){
  if(!is.null(anisotropy)){
    anisotropy = try(as.AnisotropyScaling(anisotropy))
    if(class(anisotropy)=="try-error") stop("krigingNeighbourhood: anisotropy description provided cannot be parsed")
  }
  # here space for checking that parameters are sensible
  
  # return output
  res = list(nmax=nmax, nmin=nmin, omax=omax, maxdist=maxdist, force=force, anisotropy=anisotropy)
  class(res) = "gmKrigingNeighbourhood"
  return(res)
}





##### Direct sampling parameters ---------------
#' Create a parameter set specifying a direct sampling algorithm
#'
#' Create a parameter set describing a direct sampling algorithm to multipoint simulation.
#' All parameters except `nsim` are optional, as they have default values reasonable 
#' according to experience.
#'
#' @param nsim number of realisations desired (attention: current algorithm is slow, start with small values!)
#' @param scanFraction maximum fraction of the training image to be scanned on each iteration
#' @param patternSize number of observations used for conditioning the simulation
#' @param gof maximum acceptance discrepance between a data event in the training image and the conditioning data event
#' @param ... further parameters, not used
#'
#' @return an S3-list of class "gmDirectSamplingParameters" containing the six elements given as arguments 
#' to the function. This is just a compact way to provide further functions such as [predict.gmSpatialModel()]
#' with appropriate triggers for choosing a prediction method or another, in this case for triggering 
#' direct sampling.
#' @export
#' @aliases DirectSamplingParameters DSpars
#'
#' @examples
#' (dsp = DSpars(nsim=100, scanFraction=75, patternSize=6, gof=0.05))
#' ## then run predict(..., pars=dsp)
DSpars <- DirectSamplingParameters <- function(nsim=1, scanFraction=0.25, patternSize=10, gof=0.05, ...){
  ll = list(nsim=nsim, scanFraction=scanFraction, patternSize=patternSize, gof=gof, ...)
  class(ll) = "gmDirectSamplingParameters"
  return(ll)
}




#### simulation specifications --------------
#' Create a parameter set specifying a gaussian sequential simulation algorithm
#'
#' Create a parameter set describing a sequential simulation algorithm to two-point simulation,
#' mostly for covariance or variogram-based gaussian random fields.
#' 
#' @param nsim number of realisations desired
#' @param ng a neighbourhood specification, as obtained with function [KrigingNeighbourhood()]
#' @param rank currently ignored (future functionality: obtain a reduced-rank simulation)
#' @param debug.level degree of verbosity of results; negative values produce a progress bar; values can be
#' extracted from [gstat::predict.gstat()] 
#' @param ... further parameters, currently ignored
#'
#' @return an S3-list of class "gmSequentialSimulation" containing the four elements given as arguments 
#' to the function. This is just a compact way to provide further functions such as [predict.gmSpatialModel()]
#' with appropriate triggers for choosing a prediction method or another, in this case for triggering 
#' sequential Gaussian simulation.
#' @export
#'
#' @examples
#' data("jura", package="gstat")
#' X = jura.pred[,1:2]
#' summary(X)
#' Zc = jura.pred[,7:10]
#' ng_local = KrigingNeighbourhood(maxdist=1, nmin=4, omax=5, force=TRUE)
#' (sgs_local = SequentialSimulation(nsim=100, ng=ng_local, debug.level=-1))
#' ## then run predict(..., pars=sgs_local)
SequentialSimulation = function(nsim=1, ng=NULL, rank=Inf, debug.level=1, ...){
  if(is.null(ng)) warning("SequentialSimulation: local neighbourhood is required; calculations will be stopped if the spatial model object does not include it")
  res = list(nsim=nsim, ng=ng, rank=rank, debug.level=debug.level, ...)
  class(res) = "gmSequentialSimulation"
  return(res)
}



#' Create a parameter set specifying a turning bands simulation algorithm
#'
#' Create a parameter set describing a turning bands algorithm to two-point simulation,
#' mostly for covariance or variogram-based gaussian random fields.
#' 
#' @param nsim number of realisations desired
#' @param nBands number of bands desired for the decomposition of the 2D or 3D space in individual signals 
#' @param ... further parameters, currently ignored
#'
#' @return an S3-list of class "gmTurningBands" containing the few elements given as arguments 
#' to the function. This is just a compact way to provide further functions such as [predict.gmSpatialModel()]
#' with appropriate triggers for choosing a prediction method or another, in this case for triggering 
#' turning bands simulation.
#' @export
#'
#' @examples
#' (tbs_local = TurningBands(nsim=100, nBands=300))
#' ## then run predict(..., pars=tbs_local)
TurningBands = function(nsim=1, nBands=1000,  ...){
  res = list(nsim=nsim, nBands=nBands, ...)
  class(res) = "gmTurningBands"
  return(res)
}


#' Create a parameter set specifying a LU decomposition simulation algorithm
#'
#' Create a parameter set describing a Cholesky (or LU) decomposition algorithm to two-point simulation,
#' mostly for covariance or variogram-based gaussian random fields.
#' 
#' @param nsim number of realisations desired
#' @param ... further parameters, currently ignored
#'
#' @return an S3-list of class "gmCholeskyDecomposition" containing the few elements given as arguments 
#' to the function. This is just a compact way to provide further functions such as [predict.gmSpatialModel()]
#' with appropriate triggers for choosing a prediction method or another, in this case for triggering 
#' LU or Cholesky decomposition simulation.
#' @export
#'
#' @examples
#' (chols_local = CholeskyDecomposition(nsim=100, nBands=300))
#' ## then run predict(..., pars=chols_local)
CholeskyDecomposition = function(nsim=1, ...){
  res = list(nsim=nsim,  ...)
  class(res) = "gmCholeskyDecomposition"
  return(res)
}




