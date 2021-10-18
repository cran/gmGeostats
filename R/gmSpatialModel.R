#### gmSpatialModel ------
#' Conditional spatial model data container
#' 
#' This class is devised to contain a conditional spatial model, with: some conditioning data 
#' (a [sp::SpatialPointsDataFrame()]), an unconditional geospatial model (a structure with e.g. 
#' a training image; or the information defining a Gaussian random field); and eventually some
#' extra method parameters. The class extends [sp::SpatialPointsDataFrame()] and has therefore its slots,
#' plus `model` (for the unconditional model) and `parameters` (for the extra method information)
#'
#' @slot data a data.frame (or class extending it) containing the conditional data
#' @slot coords a matrix or dataframe of 2-3 columns containing the sampling locations of the conditional data
#' @slot coords.nrs see [sp::SpatialPointsDataFrame()]
#' @slot bbox see [sp::SpatialPointsDataFrame()]
#' @slot proj4string see [sp::SpatialPointsDataFrame()]
#' @slot model gmUnconditionalSpatialModel. Some unconditional geospatial model. It can be NULL. 
#' @slot parameters gmSpatialMethodParameters. Some method parameters. It can be NULL
#' @family gmSpatialModel
#'
#' @return You will seldom create the spatial model directly. Use instead the creators `make.gm*` linked below
#' @export
#' @include abstractClasses.R
#' @importClassesFrom sp SpatialPointsDataFrame
#' 
#' @examples
#' data("jura", package="gstat")
#' library(sp)
#' X = jura.pred[,1:2]
#' Zc = jura.pred[,7:13]
#' spdf = sp::SpatialPointsDataFrame(coords=X, data=Zc)
#' new("gmSpatialModel", spdf)
#' make.gmCompositionalGaussianSpatialModel(data=Zc, coords=X, V="alr")
setClass("gmSpatialModel", contains="SpatialPointsDataFrame",
         slots = list(model="gmUnconditionalSpatialModel", 
                      parameters="gmSpatialMethodParameters")
)


### creators --------
# make.gmMultilayerSpatialModel
# make.gmMultivariateSpatialModel
# make.gmAmountsSpatialModel
# make.gmUnivariateSpatialModel
# ...

#' Construct a Gaussian gmSpatialModel for regionalized multivariate data
#' 
#' Construct a regionalized multivariate data container to be used for Gaussian-based geostatistics: variogram modelling, cokriging and simulation. 
#' @param data either a data set of any data.frame similar class, or else a [sp::SpatialPointsDataFrame()] containing it
#' @param coords the coordinates of the sampling locations, if no SpatialPointsDataFrame was provided
#' @param model a variogram model, of any relevant class
#' @param beta (see `formula`) the coefficients of dependence of the mean of the random field, if these are known; e.g. if `formula=~1` constant mean,
#' and the mean is indeed known, `beta` would be a compositional mean; seldom used directly
#' @param formula a formula without left-hand-side term, e.g. ~1 or ~Easting+Northing, specifying what do we know of the
#' dependence of the mean of the random field; this follows the same ideas than in [gstat::gstat()]
#' @param ng optional neighborhood information, typically created with [KrigingNeighbourhood()]
#' @param nmax optional, neighborhood description: maximum number of data points per cokriging system
#' @param nmin optional, neighborhood description: minimum number of data points per cokriging system
#' @param omax optional, neighborhood description: maximum number of data points per cokriging system per quadrant/octant
#' @param maxdist optional, neighborhood description: maximum radius of the search neighborhood
#' @param force optional logical, neighborhood description: if not `nmin` points are found inside `maxdist` radius, 
#' keep searching. This and all preceding arguments for neighborhood definition are borrowed from [gstat::gstat()]
#'
#' @return A "gmSpatialModel" object with all information provided appropriately structured. See [gmSpatialModel-class].
#' @export
#' @family gmSpatialModel 
#' @seealso [SequentialSimulation()], [TurningBands()] or [CholeskyDecomposition()] for specifying the exact 
#' simulation method and its parameters, [predict_gmSpatialModel] for running predictions or simulations
#'
#' @examples
#' data("jura", package="gstat")
#' X = jura.pred[,1:2]
#' Zc = jura.pred[,7:13]
#' make.gmMultivariateGaussianSpatialModel(data=Zc, coords=X)
make.gmMultivariateGaussianSpatialModel <- function(
  data, coords = attr(data, "coords"), ## data trench
  model=NULL, beta=model$beta, formula=model$formula,  ## model trench
  ng=NULL, nmax=ng$nmax, nmin=ng$nmin, omax=ng$omax, maxdist=ng$maxdist, force=ng$force){ # method trench
  # consider the class of 'data' 
  if(is(data,"Spatial")){
    dt = data@data
  }else{
    dt = data
  }
  spdf = sp::SpatialPointsDataFrame(sp::SpatialPoints(coords), data=dt)
  ng = KrigingNeighbourhood(nmax=c(nmax,Inf)[1], nmin=c(nmin,0)[1], omax=c(omax,0)[1], maxdist=c(maxdist,Inf)[1], force=c(force,FALSE)[1])
  if(is.null(formula)) formula=~1
  if(is.null(beta)) beta=Inf
  if(is.null(model)){
    vgm = new("gmGaussianModel", structure=NULL, beta=beta, formula=formula)
  }else{
    vgm = new("gmGaussianModel", structure=model, beta=beta, formula=formula)
  }
  res = new("gmSpatialModel", spdf, model=vgm, parameters=ng)
  return(res)
}


#' Construct a Gaussian gmSpatialModel for regionalized compositions
#' 
#' Construct a regionalized compositional data container to be used for Gaussian-based geostatistics: variogram modelling, cokriging and simulation. 
#' @param data either a [compositions::acomp()] compositional data set, or else a [sp::SpatialPointsDataFrame()] containing it
#' @param coords the coordinates of the sampling locations, if no SpatialPointsDataFrame was provided
#' @param V optionally, a matrix of logcontrasts, or else one of the following strings: "alr", "ilr" or "clr"; 
#' to produce a plot of the empirical variogram in the corresponding representation; default to variation-variograms
#' @param prefix the desired prefix name for the logratio variables, if this is wished to be forced; otherwise derived from `V`
#' @param model a variogram model, of any relevant class
#' @param beta (see `formula`) the coefficients of dependence of the mean of the random field, if these are known; e.g. if `formula=~1` constant mean,
#' and the mean is indeed known, `beta` would be a compositional mean; seldom used directly
#' @param formula a formula without left-hand-side term, e.g. `~1` or `~Easting+Northing`, specifying what do we know of the
#' dependence of the mean of the random field; this follows the same ideas than in [gstat::gstat()]
#' @param ng optional neighborhood information, typically created with [KrigingNeighbourhood()]
#' @param nmax optional, neighborhood description: maximum number of data points per cokriging system
#' @param nmin optional, neighborhood description: minimum number of data points per cokriging system
#' @param omax optional, neighborhood description: maximum number of data points per cokriging system per quadrant/octant
#' @param maxdist optional, neighborhood description: maximum radius of the search neighborhood
#' @param force optional logical, neighborhood description: if not `nmin` points are found inside `maxdist` radius, 
#' keep searching. This and all preceding arguments for neighborhood definition are borrowed from [gstat::gstat()]
#'
#' @return A "gmSpatialModel" object with all information provided appropriately structured. See [gmSpatialModel-class].
#' @export
#' @family gmSpatialModel
#' @seealso [SequentialSimulation()], [TurningBands()] or [CholeskyDecomposition()] for specifying the exact 
#' simulation method and its parameters, [predict_gmSpatialModel] for running predictions or simulations
#'
#' @examples
#' data("jura", package="gstat")
#' X = jura.pred[1:20,1:2]
#' Zc = compositions::acomp(jura.pred[1:20,7:13])
#' make.gmCompositionalGaussianSpatialModel(data=Zc, coords=X, V="alr")
make.gmCompositionalGaussianSpatialModel <- function(
  data, coords = attr(data, "coords"), V="ilr", prefix=NULL,  ## data trench
  model=NULL, beta=model$beta, formula=model$formula,  ## model trench
  ng=NULL, nmax=ng$nmax, nmin=ng$nmin, omax=ng$omax, maxdist=ng$maxdist, force=ng$force)  # method trench
{
  if(is(data,"Spatial")){
    dt = data@data
  }else{
    dt = data
  }
  if(is.rmult(dt)){
    W = gsi.getV(dt)
    V = t(gsiInv(W))
    warning("make.gmCompositionalSpatialModel: data of class 'rmult'! provided V will be ignored")
    spdf = sp::SpatialPointsDataFrame(sp::SpatialPoints(coords), data=dt)
  }else{
    o = gsi.produceV(V=V,D=ncol(dt),orignames=colnames(dt),giveInv=FALSE, prefix=prefix)
    prefix = o$prefix
    V = o$V
    if(is.null(rownames(V))) tryCatch(rownames(V) <- colnames(data))
    if(is.null(rownames(V))) rownames(V) = paste("v", 1:nrow(V), sep="")
    if(is.null(colnames(V))) colnames(V) = paste(prefix, 1:ncol(V), sep="")
    spdf = sp::SpatialPointsDataFrame(sp::SpatialPoints(coords), data=ilr(dt,V))
  }
  if(is.null(ng)) ng = KrigingNeighbourhood(nmax=c(nmax,Inf)[1], nmin=c(nmin,0)[1], omax=c(omax,0)[1], maxdist=c(maxdist,Inf)[1], force=c(force,FALSE)[1])
  if(is.null(formula)) formula=~1
  if(is.null(beta)) beta=Inf
  if(is.null(model)){
    vgm = new("gmGaussianModel", structure=NULL, beta=beta, formula=formula)
  }else{
    vgm = new("gmGaussianModel", structure=model, beta=beta, formula=formula)
  }
  res = new("gmSpatialModel", spdf, model=vgm, parameters=ng)
  return(res)
} 





#' Construct a Multi-Point gmSpatialModel for regionalized compositions
#' 
#' Construct a regionalized compositional data container to be used for multipoint geostatistics. 
#' @param data either a [compositions::acomp()] compositional data set, or else a [sp::SpatialPointsDataFrame()] containing it
#' @param coords the coordinates of the sampling locations, if no SpatialPointsDataFrame was provided
#' @param V optionally, a matrix of logcontrasts, or else one of the following strings: "alr", "ilr" or "clr"; 
#' to produce a plot of the empirical variogram in the corresponding representation; default to variation-variograms
#' @param prefix the desired prefix name for the logratio variables, if this is wished to be forced; otherwise derived from `V`
#' @param model a training image, of any appropriate class (typically a [sp::SpatialGridDataFrame()] or [sp::SpatialPixelsDataFrame()]) 
#'
#' @return  A "gmSpatialModel" object with all information provided appropriately structured. See [gmSpatialModel-class].
#' @export
#' @family gmSpatialModel 
#' @seealso [DirectSamplingParameters()] for specifying a direct simulation method parameters,
#' [predict_gmSpatialModel] for running the simulation
make.gmCompositionalMPSSpatialModel = function(
  data, coords = attr(data, "coords"), V="ilr", prefix=NULL,  ## data trench
  model=NULL)   ## model trench
{
  # extract data
  if(is(data,"Spatial")){
    dt.dt = data@data
  }else{
    dt.dt = data
  }
  # extract grid topology
  if("grid" %in% slotNames(data)){
    grid.dt = data@grid
  }else if("grid" %in% slotNames(coords)){
    grid.dt = coords@grid
  }else{
    grid.dt = NULL
    warning("make.gmCompositionalMPSSpatialModel: grid topology is missing, it will be inferred from the data")
  } 
  # check the adequacy of the model provided
  if(is(model,"Spatial")){
    dt.ti = model@data
    grid.ti = if("grid" %in% slotNames(model)) model@grid else NULL
    coords.ti = if(is(model, "SpatialGridDataFrame")) NULL else sp::coordinates(model)
  }else if(is.data.frame(model) | is.array(model)){
    warning("make.gmCompositionalMPSSpatialModel: model provided is not a Spatial object; conversion will be attempted")
    dt.ti = model
    grid.ti = NULL
    coords.ti = attr(model, "coords")
  }else if(is.null(model)){
    dt.ti = dt.dt
    grid.ti = NULL
    coords.ti = NULL
  }else stop("make.gmCompositionalMPSSpatialModel: model provided cannot be interpreted")
  # create grid for the TI or else check consistency
  if(!is.null(grid.dt)){
    if(is.null(grid.ti)){
      x0 = grid.dt@grid@cellcentre.offset
      Dx = grid.dt@grid@cellsize
      n = dim(dt.ti)[1:2]
      grid.ti = sp::GridTopology(cellcentre.offset = x0, cellsize = Dx, cells.dim = n)
    }else{
      stopifnot(all(grid.dt@cellsize==grid.ti@cellsize))
    }
  }else{
    if(!is.null(grid.ti)){
      x0 = grid.ti@grid@cellcentre.offset
      Dx = grid.ti@grid@cellsize
      n = grid.ti@grid@cellsize
      grid.dt = sp::GridTopology(cellcentre.offset = x0, cellsize = Dx, cells.dim = n)
      warning("make.gmCompositionalMPSSpatialModel: grid topology found on the TI but not on the data; TI from the grid will be used for the data!")
    }
  }
  # check consistency between data parts of the data and model
  if(ncol(dt.dt)!=ncol(dt.ti)) stop("make.gmCompositionalMPSSpatialModel: number of variables for the data and model provided do not coincide!")
  # compute object data
  if(is.rmult(dt.dt)){
    W = gsi.getV(dt.dt)
    V = t(gsiInv(W))
    warning("make.gmCompositionalMPSSpatialModel: data of class 'rmult'! provided V will be ignored")
    spdf = sp::SpatialPixelsDataFrame(sp::SpatialPoints(coords), data=dt.dt, grid=grid.dt)
  }else{
    o = gsi.produceV(V=V,D=ncol(dt.dt),orignames=colnames(dt.dt),giveInv=FALSE, prefix=prefix)
    prefix = o$prefix
    V = o$V
    if(is.null(rownames(V))) tryCatch(rownames(V) <- colnames(dt.dt))
    if(is.null(rownames(V))) rownames(V) = paste("v", 1:nrow(V), sep="")
    if(is.null(colnames(V))) colnames(V) = paste(prefix, 1:ncol(V), sep="")
    spdf = sp::SpatialPixelsDataFrame(sp::SpatialPoints(coords), 
                                      data=compositions::ilr(dt.dt,V), grid=grid.dt)
  }
  myilr = function(x,V){
    res = log(as.matrix(x)) %*% V
    tk = gmApply(is.na(x),1, any)
    res[tk,] =NA
    return(data.frame(res))
  } 
  # compute object model
  if(is(model, "SpatialGridDataFrame")){
    model = sp::SpatialGridDataFrame(grid = sp::getGridTopology(model), 
                                   data = myilr(model@data,V))
  }else if(is(model, "SpatialPixelsDataFrame")){
    model = sp::SpatialPixelsDataFrame(points = 
                                         sp::SpatialPoints(sp::coordinates(model)), 
                                   grid = sp::getGridTopology(model), 
                                   data = myilr(model@data,V))
  }else if(is.null(model)){
  }else if(is.null(coords.ti)){
    
    model = sp::SpatialGridDataFrame(grid = grid.ti, data = myilr(dt.ti,V))
  }else{
    model = sp::SpatialPixelsDataFrame(points = sp::SpatialPoints(coords.ti), 
                                       grid = grid.ti, data = myilr(dt.ti,V))
  }
  res = new("gmSpatialModel", spdf, model=model)
  return(res)
}


### exporter to gstat -----
as.gstat.gmSpatialModel <- function(object, ...){
  # extra arguments
  lldots = list(...)
  # data elements
  coords = sp::coordinates(object)
  X = compositions::rmult(object@data, V= gsi.getV(object@data), orig=gsi.orig(object@data))
  V = gsi.getV(X)
  if(!is.null(V)){
    # compositional case
    compo = backtransform(X)
    Vinv = t(gsiInv(V))
    prefix = sub("1","",colnames(V)[1])
    if(is.null(prefix) | length(prefix)==0)   prefix = sub("1","",colnames(object@data)[1])
    if(is.null(prefix) | length(prefix)==0)   prefix = "ilr"
    # model elements
    if(!is(object@model, "gmGaussianModel")) stop("as.gstat: object@model must be of class 'gmGaussianModel'!")
    if(is.null(aux <- object@model@structure)){
      lrvgLMC = NULL
    }else{
      lrvgLMC = as.LMCAnisCompo(aux, V=V) 
    }
    formulaterm = paste(as.character(object@model@formula), collapse="")
    beta = object@model@beta
    if(any(is.infinite(beta))) beta = NULL
    # neighbourhood
    ng = object@parameters
    # manage manual changes of parameters given in dots...
    for(nm in names(ng)){
      tk = grep(nm, names(lldots))
      if(length(tk)>1) ng[[tk]] = lldots[[tk]] 
    }
    # convert
    if(!is(ng, "gmKrigingNeighbourhood")) stop("as.gstat: object@parameters must be of class 'gmGaussianMethodParameters'!")
    res = compo2gstatLR(coords=coords, compo=compo, V=Vinv, lrvgLMC=lrvgLMC, 
                        nscore=FALSE, formulaterm = formulaterm, prefix=prefix, beta=beta, 
                        nmax=ng$nmax, nmin=ng$nmin, omax=ng$omax, maxdist=ng$maxdist, force=ng$force, 
                        ...)
    return(res)
  }else{
    # non-compositional case
    prefix = sub("1","",colnames(V)[1])
    vgLMC <- object@model@structure
    formulaterm = paste(as.character(object@model@formula), collapse="")
    beta = object@model@beta
    if(any(is.infinite(beta))) beta = NULL
    # neighbourhood
    ng = object@parameters
    # manage manual changes of parameters given in dots...
    for(nm in names(ng)){
      tk = grep(nm, names(lldots))
      if(length(tk)>1) ng[[tk]] = lldots[[tk]] 
    }
    res = rmult2gstat(coords=coords, data=X, V=V, vgLMC=vgLMC, 
                    nscore=FALSE, formulaterm = formulaterm, prefix=prefix, beta=beta, 
                    nmax=ng$nmax, nmin=ng$nmin, omax=ng$omax, maxdist=ng$maxdist, force=ng$force,
                    ...)    
  }
}



#' Recast spatial object to gmSpatialModel format
#' 
#' Recast a spatial data object model to format gmSpatialModel
#'
#' @param object object to recast
#' @param ... extra parameters for generic functionality
#'
#' @return The same spatial object re-structured as a "gmSpatialModel", see [gmSpatialModel-class]
#' @export
#' @family gmSpatialModel
as.gmSpatialModel <- function(object, ...) UseMethod("as.gmSpatialModel", object)

#' @describeIn as.gmSpatialModel Recast spatial object to gmSpatialModel format
#' @method as.gmSpatialModel default
as.gmSpatialModel.default = function(object, ...) object


#' @describeIn as.gmSpatialModel Recast spatial object to gmSpatialModel format
#' @param V optional, if the original data in the sptail object was compositional, which logcontrasts 
#' were used to express it? Thsi can be either one string of "alr", "ilr" or "clr", or else a 
#' (Dx(D-1))-element matrix of logcontrasts to pass from compositions to logratios
#' @method as.gmSpatialModel gstat
as.gmSpatialModel.gstat = function(object, V=NULL, ...){
  stop("as.gmSpatialModel.gstat: not yet implemented")
}


#' Predict method for objects of class 'gmSpatialModel'
#' 
#' This is a one-entry function for several spatial prediction and simulation methods, for model objects 
#' of class [gmSpatialModel-class]. The several methods are chosen by means of `pars` objects of the 
#' appropriate class. 
#'
#' @param object a complete "gmSpatialModel", containing conditioning data and unconditional model 
#' @param newdata a collection of locations where a prediction/simulation is desired; this is typically 
#' a [sp::SpatialPoints()], a data.frame or similar of X-Y(-Z) coordinates; or perhaps for gridded data 
#' an object of class  [sp::GridTopology()], [sp::SpatialGrid()] or [sp::SpatialPixels()]
#' @param pars parameters describing the method to use, *enclosed in an object of appropriate class* 
#' (according to the method; see below)
#' @param ... further parameters for generic functionality, currently ignored
#'
#' @return Depending on the nature of `newdata`, the result will be a data container of the same kind, 
#' extended with the predictions or simulations. For instance, if we want to obtain predictions on the 
#' locations of a "SpatialPoints", the result will be a [sp::SpatialPointsDataFrame()]; if we want to obtain
#' simulations on the coordinates provided by a "data.frame", the result will be a [DataFrameStack()] with 
#' the spatial coordinates stored as an extra attribute; or if the input for a simulation is a masked grid of class
#' [sp::SpatialPixels()], the result will be of class [sp::SpatialPixelsDataFrame()] which `data` slot will be 
#' a [DataFrameStack]. 
#' 
#' @details Package "gmGeostats" aims at providing a broad series of algorithms for geostatistical prediction
#' and simulation. All can be accesses through this interface, provided that arguments `object` and `pars` are of the 
#' appropriate kind. In `object`, the most important criterion is the nature of its slot `model`. In `pars`
#' its class counts: for the creation of informative parameters in the appropriate format and class, a series
#' of accessory functions are provided as well.
#' 
#' Classical (gaussian-based two-point) geostatistics are obtained if `object@model` contains a covariance function,
#' or a variogram model. Argument `pars` can be created with functions such as [KrigingNeighbourhood()],
#' [SequentialSimulation()], [TurningBands()] or [CholeskyDecomposition()] to respectively trigger a cokriging, as
#' sequential Gaussian simulation, a turning bands simulation, or a simulation via Cholesky decomposition.
#' The kriging neighbourhood can as well be incorporated in the "gmSpatialModel" `object` directly, or even be
#' nested in a "SequentialSimulation" parameter object.
#' 
#' Conversely, to run a multipoint geostatistics algorithm, the first condition is that `object@model` contains a 
#' training image. Additionally, `pars` must describe the characteristics of the algorithm to use. Currently, only
#' direct sampling is available: it can be obtained by providing some parameter object created with a call to
#' [DirectSamplingParameters()]. This method requires `newdata` to be on a gridded set of locations (acceptable
#' object classes are `sp::gridTopology`, `sp::SpatialGrid`, `sp::SpatialPixels`, `sp::SpatialPoints` or `data.frame`,
#' for the last two a forced conversion to a grid will be attempted).  
#' @family gmSpatialModel
#' @name predict_gmSpatialModel
NULL



#' @rdname predict_gmSpatialModel
#' @export
#' @method predict gmSpatialModel
predict.gmSpatialModel <- function(object, newdata, pars=NULL, ...){
  if(is.null(pars)){
    return(Predict(object, newdata, ...))
  }else{
    return(Predict(object, newdata, pars, ...))
  }
}

#' @rdname predict_gmSpatialModel
#' @export
setMethod("predict", signature(object="gmSpatialModel"), definition = predict.gmSpatialModel)



#' @rdname predict_gmSpatialModel
#' @include gmAnisotropy.R
#' @include preparations.R
#' @export
setMethod("Predict",signature(object="gmSpatialModel", newdata="ANY"),
          function(object, newdata, ...){
            if(is.null(object@parameters)) object@parameters = KrigingNeighbourhood() 
            Predict(object, newdata, pars = object@parameters , ...)
          }
)



#' @rdname predict_gmSpatialModel
#' @include gmSpatialMethodParameters.R
#' @export
setMethod("Predict",signature(object="gmSpatialModel", newdata="ANY", pars="gmNeighbourhoodSpecification"),
          function(object, newdata, pars, ...){
            cat("starting cokriging \n")
            object@parameters = pars
            out = predict(as.gstat(object), newdata=newdata, ...)
            return(out)
          }
)



#' @rdname predict_gmSpatialModel
#' @export
setMethod("Predict",signature(object="gmSpatialModel", newdata="ANY", pars="gmTurningBands"),
          function(object, newdata, pars, ...){
            stop("Turning Bands method not yet interfaced here; use")
            cat("starting turning bands \n")
          }
)


#' @rdname predict_gmSpatialModel
#' @export
setMethod("Predict",signature(object="gmSpatialModel", newdata="ANY", pars="gmCholeskyDecomposition"),
          function(object, newdata, pars, ...){
            stop("Choleski decomposition method not yet implemented")
            cat("starting Choleski decomposition \n")
          }
)




#' @rdname predict_gmSpatialModel
#' @export
setMethod("Predict",signature(object="gmSpatialModel", newdata="ANY", pars="gmSequentialSimulation"),
          function(object, newdata, pars, ...){
            cat("starting SGs \n")
            object@parameters = pars$ng
            erg = predict(as.gstat(object), newdata=newdata, nsim=pars$nsim, debug.level=pars$debug.level, ...)
            Dg = ncol(object@coords)
            erg = DataFrameStack(erg[,-(1:Dg)], 
                                 dimnames=list(
                                   loc=1:nrow(erg), sim=1:pars$nsim,
                                   var=colnames(object@data)
                                 ),
                                 stackDim="sim")
            attr(erg, "coords") = newdata
            return(erg)
          }
)



#' @rdname predict_gmSpatialModel
#' @export
setMethod("Predict",signature(object="gmSpatialModel", newdata="ANY", pars="gmDirectSamplingParameters"),
          function(object, newdata, pars, ...){
            
            cat("starting direct sampling \n")
            # extract training image
            gt.ti = sp::getGridTopology(object@model)
            dt.ti = as(object@model,"SpatialGridDataFrame")@data
            # extract newdata mask
            mask = getMask(newdata)
            # fuse newdata, conditioning data and mask
            if(is.null(newdata)){
              gt.nw = NULL # object@grid
              newdata = sp::SpatialPixelsDataFrame(points = as(object, "SpatialPoints"), data=object@data)  
              # tolerance = sqrt(sum(gt.nw@cellsize^2))/2 , grid=gt.nw)
            }else if(is(newdata, "GridTopology")){
              gt.nw = newdata
              newdata = sp::SpatialPixelsDataFrame(points = as(object, "SpatialPoints"), data=object@data, 
                                                   tolerance = sqrt(sum(gt.nw@cellsize^2))/2 , grid=gt.nw)
            }else if(is(newdata, "SpatialGrid")){
              gt.nw = sp::getGridTopology(newdata)
              newdata = sp::SpatialPixelsDataFrame(points = as(object, "SpatialPoints"), data=object@data, 
                                                   tolerance = sqrt(sum(gt.nw@cellsize^2))/2 , grid=gt.nw)
            }else if(is(newdata, "SpatialPixels")){
              gt.nw = sp::getGridTopology(newdata)
              cc = rbind(sp::coordinates(newdata), sp::coordinates(object))
              ### WARNING: This is how it works!! COrrect from SpatialPoints and data.frame
              aux = matrix(NA, nrow=nrow(sp::coordinates(newdata)), ncol=ncol(object@data))
              colnames(aux) = colnames(object@data)
              dt = rbind( aux , object@data)
              newdata = sp::SpatialPixelsDataFrame(points = sp::SpatialPoints(cc), data=dt, 
                                                   tolerance = sqrt(sum(gt.nw@cellsize^2))/2 , grid=gt.nw)
            }else if(is(newdata, "SpatialPoints")){
              gt.nw = object@model@grid
              cc = rbind(sp::coordinates(object), sp::coordinates(newdata))
              ### WARNING: NOT YET TESTED WHAT HAPPENS WITH DUPLICATE LOCATIONS!!! 
              aux = matrix(NA, nrow=nrow(sp::coordinates(newdata)), ncol=ncol(object@data))
              colnames(aux) = colnames(object@data)
              dt = rbind( object@data, aux )
              newdata = sp::SpatialPixelsDataFrame(points = sp::SpatialPoints(cc), data=dt,
                                                   tolerance = sqrt(sum(gt.nw@cellsize^2))/2) # , grid=gt.nw)
            }else if(is.data.frame(newdata)){
              gt.nw = NULL # object@grid
              cc = rbind(sp::coordinates(object), newdata[,1:2])
              if( all( colnames(object@data) %in% colnames(newdata) ) ){
                aux = as.matrix(newdata[,colnames(object@data)])
              }else{
                aux = matrix(NA, nrow=nrow(sp::coordinates(newdata)), ncol=ncol(object@data))
              }
              colnames(aux) = colnames(object@data)
              ### WARNING: NOT YET TESTED WHAT HAPPENS WITH DUPLICATE LOCATIONS!!! 
              dt = rbind( object@data, aux )
              newdata = sp::SpatialPixelsDataFrame(points = sp::SpatialPoints(cc), data=dt)
              # tolerance = sqrt(sum(gt.nw@cellsize^2))/2 , grid=gt.nw)
            }
            if(is.null(gt.nw)) gt.nw = sp::getGridTopology(newdata)
            dt.nw = as(newdata,"SpatialGridDataFrame")@data
            dt.ti = as(object@model,"SpatialGridDataFrame")@data
            if(all(gt.nw@cellsize!=gt.ti@cellsize)) 
              stop("predict for gmSpatialModel with gmDirectSamplingParameters: inferred grid topologies for newdata, conditioning data and model do not coincide")
            if(is.null(mask)) mask = rep(TRUE, nrow(dt.nw))
            #erg = gsi.DS4CoDa(n=pars$patternSize, f=pars$scanFraction, t=pars$gof, n_realiz=pars$nsim, 
            #            nx_TI=gt.ti@cells.dim[1], ny_TI=gt.ti@cells.dim[2], 
            #            nx_SimGrid= gt.nw@cells.dim[1], ny_SimGrid=gt.nw@cells.dim[2],
            #            TI_input=as.matrix(dt.ti),
            #            SimGrid_input=as.matrix(dt.nw), 
            #            V = "I", W=gsi.getV(object@data), 
            #            ivars_TI = colnames(dt.ti), 
            #            SimGrid_mask = mask, 
            #            invertMask = FALSE)
            erg = gsi.DS(n=pars$patternSize, f=pars$scanFraction, t=pars$gof, n_realiz=pars$nsim, 
                         dim_TI=gt.ti@cells.dim, dim_SimGrid=gt.nw@cells.dim,   
                         TI_input=as.matrix(dt.ti), SimGrid_input=as.matrix(dt.nw),
                         ivars_TI = colnames(dt.ti), 
                         SimGrid_mask = mask, 
                         invertMask = FALSE)
            erg = gmApply(erg, FUN=ilrInv, V=gsi.getV(object@data), orig=gsi.orig(object@data))
            erg = sp::SpatialGridDataFrame(grid = gt.nw, data=erg)
            return(erg)
          }
)


#' @describeIn gmSpatialModel Compute a variogram, see [variogram_gmSpatialModel()] for details
#' @inheritParams variogram_gmSpatialModel
#' @include variograms.R
#' @export
setMethod("variogram", signature=c(object="gmSpatialModel"), 
          def=variogram_gmSpatialModel)



#' @describeIn gmSpatialModel S4 wrapper method around [logratioVariogram()] for `gmSpatialModel`
#' objects
#' @inheritParams logratioVariogram
#' @inheritParams logratioVariogram_gmSpatialModel
#' @include compositionsCompatibility.R
#' @export
setMethod("logratioVariogram", "gmSpatialModel", logratioVariogram_gmSpatialModel)


#' @describeIn gmSpatialModel convert from gmSpatialModel to gstat; see [as.gstat()]
#' for details
#' @inheritParams as.gstat
#' @export
#' @include gstatCompatibility.R
setMethod("as.gstat", signature="gmSpatialModel", def=as.gstat.gmSpatialModel)



### tests -----------
# if(!exists("do.test")) do.test=FALSE
# if(do.test){
#   library("gstat")
#   library("magrittr")
#   data("jura", package = "gstat")
#   dt = jura.pred %>% dplyr::select(Cd:Zn)
#   X = jura.pred[,1:2]
#   a1 = make.gmCompositionalGaussianSpatialModel(acomp(dt), X)
#   spc = SpatialPointsDataFrame(SpatialPoints(X),acomp(dt))
#   a2 = make.gmCompositionalGaussianSpatialModel(spc)
#   a3 = make.gmCompositionalGaussianSpatialModel(spc, V="alr")
#   a3gs = as.gstat(a3)
#   a3vg = variogram(a3gs) 
#   plot(a3vg)
#   a3gs = fit_lmc(a3vg, a3gs, vgm("Sph", psill=1, nugget=0, range=1))
#   plot(a3vg, a3gs$model)
#   a4 = make.gmCompositionalGaussianSpatialModel(spc, V="alr", model=a3gs$model)
#   a4gs = as.gstat(a4)
#   plot(a3vg, a4gs$model)
# }





### converters -----
# gmSpatialModel2SpatialGridDataFrame = function(from, to){
#   tol = ifelse(is(grid, "GridTopology"), 0.5*sqrt(sum(from@grid@cellsize^2)), sqrt(.Machine$double.eps))
#   to = SpatialPixelsDataFrame(
#     grid = from@grid, data = from@data, points = SpatialPoints(from@coords), 
#     proj4string = proj4string(from), tolerance = tol
#   )
#   return(as(to, "SpatialGridDataFrame"))
# }
# SpatialGridDataFrame4gmSpatialModel = function(from, value){
#   gt.mdl = from@model@grid
#   gt.new = value@grid
#   if(!is.null(gt.mdl))
#     if(!all(gt.new@cellsize==gt.mdl@cellsize)) stop("replace(SpatialGridDataFrame->gmSpatialModel): provided sgdf grid topology inconsistent with existing model topology")
#   from@grid = gt.new
#   from@data = value@data
#   from@coords = sp::coordinates(value)
#   from@bbox = bbox(value)
#   from@proj4string = proj4string(value)
#   return(from)
# }
#setIs(class1 = "gmSpatialModel", class2 ="SpatialGridDataFrame",
#      coerce = gmSpatialModel2SpatialGridDataFrame, 
#      replace = SpatialGridDataFrame4gmSpatialModel)
#setIs(class1 = "gmSpatialModel", class2 ="SpatialPointsDataFrame",
#      coerce = function(from, to){
#        to = SpatialPointsDataFrame(
#          coords = SpatialPoints(from@coords), data = from@data,  
#          proj4string = proj4string(from), bbox = bbox(from)
#        )
#        return(to)
#      }, replace = function(from, value){
#        from@grid = NULL
#        from@data = value@data
#        from@coords = sp::coordinates(value)
#        from@bbox = bbox(value)
#        from@proj4string = proj4string(value)
#        return(from)
#      }
#)


