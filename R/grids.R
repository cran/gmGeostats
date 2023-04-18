




#' Set or get the ordering of a grid
#' 
#' Specify or retrieve the ordering in which a grid is stored in a vector (or matrix).
#'
#' @param x a data container for the elements of the grid; the grid order is stored as an attribute to it
#' @param refpoint a string specifying which point of the grid corresponds 
#' to the first element of `x`; see below
#' @param cycle a permutation of the integers `1:G` (see below)
#'
#' @return `setGridOrder(x,...)` returns the object `x` with the grid order description attached
#' as an attribute "gridOrder"; `getGridOrder(x)` retrieves this attribute and returns it.
#' @details A "gridOrder" attribute is a list consisting of two named elements:
#' \describe{
#'   \item{refpoint}{one of  "topleft", "bottomleft", "topright" or "bottomright" in 2D, or 
#'   also of "topleftsurf", "bottomleftsurf", "toprightsurf", "bottomrightsurf", "topleftdeep","bottomleftdeep",
#'   "toprightdeep" or "bottomrightdeep" in 3D ("deep" is accessory, i.e.  "topleft"== "topleftdeep"),
#'   indicating the location on the grid of the first point of the object `x`}
#'   \item{cycle}{a permutation of `1:G` indicating in which order run the dymensions, from faster to slower}
#' }
#' 
#' Thus, a conventional ordering of a (nX*nY)-element vector into a matrix to plot with [graphics::image()]
#' corresponds to an `refpoint="bottomleft"` and `cycle=1:2`, i.e. start with the lower left corner 
#' and run first by rows (eastwards), then by columns (northwards). This  is constructed by 
#' `gridOrder_array(G)`, and can be directly set to an object `x` by `setGridOrder_array(x,G)`;
#' `gridOrder_GSLib` is an alias for `gridOrder_array`.
#' 
#' The grids from package "sp" (and many other in R), on the contrary follow the convention 
#' `refpoint="topleft"` and `cycle=1:2`,  i.e. start with the upper left corner 
#' and run first by rows (eastwards), then by columns (**southwards**).  This  is constructed by 
#' `gridOrder_sp(G)`, and can be directly set to an object `x` by `setGridOrder_sp(x,G)`; `gridOrder_gstat`
#' is an alias for `gridOrder_sp`.
#' 
#' @export
#' @aliases getGridOrder
#' @seealso [sortDataInGrid()] for ways of reordering a grid
#' @examples
#' gt = sp::GridTopology(cellcentre.offset=c(1,11), cellsize=c(1,1), cells.dim=c(5,3))
#' sp::coordinates(sp::SpatialGrid(grid=gt))
#' gridOrder_sp(2)
setGridOrder = function(x, refpoint, cycle){
  refpoints =  gmApply(expand.grid(c("top","bottom"), c("left","right"), c("","surf", "deep")),1, paste, collapse="")
  if(!(refpoint %in% refpoints)) stop(paste("setGridOrder: refpoint not understood, must be one of ", paste(refpoints, collapse=", ")))
  if(any(sort(cycle)!=(1:length(cycle)))) stop("setGridOder: cycle must be a permutation of c(1,2) in 2D, or c(1,2,3) in 3D")
  attr(x,"gridOrder") = list(refpoint=refpoint, cycle=cycle)
  return(x)
}

#' @export
getGridOrder = function(x){
  res = attr(x, "gridOrder")
  if(is.null(res)){
    warning("getGridOrder: no grid ordering info found, taking the default, from options()$gmGeostats$gridOrder")
    res = options()$gmGeostats$gridOrder
  }
  return(res)
}


#' @describeIn setGridOrder Set or get the ordering of a grid
#' @param G number of geographic dimensions of the setting, typically `G=2` or `G=3`
#' @export
setGridOrder_sp = function(x, G=2){
  attr(x,"gridOrder") = gridOrder_gstat(G=G)
  return(x)
}


#' @describeIn setGridOrder Set or get the ordering of a grid
#' @export
setGridOrder_array = function(x, G=2){
  attr(x,"gridOrder") = gridOrder_array(G=G) 
  return(x)
}

#' @describeIn setGridOrder Set or get the ordering of a grid
#' @export
gridOrder_sp <-   function(G=2){
  refpoint = ifelse(G==2, "topleft", "topleftdeep")
  cycle = 1:G
  return(list(refpoint=refpoint, cycle=cycle))
}


#' @describeIn setGridOrder Set or get the ordering of a grid
#' @export
gridOrder_gstat <- gridOrder_sp
  

#' @describeIn setGridOrder Set or get the ordering of a grid
#' @export
gridOrder_array <- function(G=2){
  refpoint = ifelse(G==2, "bottomleft", "bottomleftdeep")
  cycle = 1:G
  return(list(refpoint=refpoint, cycle=cycle))
}


#' @describeIn setGridOrder Set or get the ordering of a grid
#' @export
gridOrder_GSLib <-gridOrder_array



#### internal gstatCokriging2something function ------
#' Reorganisation of cokriged compositions  
#' 
#' Produce compositional predictions out of a [gstat::gstat()] prediction  
#'
#' @param COKresult output of a [gstat::predict.gstat()] cokriging, 
#' typically of class "data.frame", [sp::SpatialPointsDataFrame()], 
#' [sp::SpatialGridDataFrame()] or [sp::SpatialPixelsDataFrame()]
#' @param ... further arguments needed for nscore (**deprecated**) 
#'
#' @return an (N,D)-object of class `c("spatialGridAcomp","acomp")` 
#' with the predictions, together with an extra attribute "krigVar" 
#' containing the cokriging covariance matrices  in an (N, D, D)-array; here N=number of 
#' interpolated locations, D=number of original components of the composition
#' @aliases gsi.gstatCokriging2rmult
#' @export
#' @seealso [image_cokriged.spatialGridRmult()] for an example
gsi.gstatCokriging2compo <- function(COKresult, ...) UseMethod("gsi.gstatCokriging2compo", COKresult)


#' @describeIn gsi.gstatCokriging2compo Reorganisation of cokriged compositions  
#' @method gsi.gstatCokriging2compo default
#' @export
gsi.gstatCokriging2compo.default <- function(COKresult, ...){
  stopifnot(is(COKresult, "Spatial"), "data" %in% slotNames(COKresult))
  coord = sp::coordinates(COKresult)
  res = gsi.gstatCokriging2compo(cbind(as.data.frame(coord), COKresult@data), ...)
  if("grid" %in% slotNames(COKresult)){
    G = ncol(coord)
    dr = 1+(getGridTopology(COKresult)@cellsize<0)
    rp = paste(c("top","bottom")[dr[1]], c("left","right")[dr[2]], sep="")
    if(G==3) rp = paste(rp, c("deep", "surf")[dr[3]], sep="")
    res = setGridOrder(res, refpoint = rp, cycle = 1:G)
  } 
  return(res)
}


#' @describeIn gsi.gstatCokriging2compo Reorganisation of cokriged compositions  
#' @method gsi.gstatCokriging2compo data.frame
#' @param V string or matrix describing which logratio was applied ("ilr", "alr", 
#' or a matrix computing the ilr corrdinates; clr is not allowed!)
#' @param orignames names of the original components (optional, but recommended)
#' @param tol for generalized inversion of the matrix (**rarely touched!**)
#' @param nscore boolean, were the data normal score-transformed? (**deprecated**)
#' @param gg in the case that normal score transformation was applied, provide the gstat object! (**deprecated**)
#' @export
gsi.gstatCokriging2compo.data.frame = function(COKresult, # output of predict.gstat
                                               V=NULL,    # string or matrix describing which logratio was applied
                                               orignames=NULL, # names of the original components (OPTIONAL)
                                               tol=1e-12, # for generalized inversion of the matrix (RARELY USED)
                                               nscore=FALSE, # were the data NS-transformed?
                                               gg=NULL, # in the case that NS was applied, provide the gstat object!
                                               ... # further arguments needed for nscore
){
  if(is.null(V)) stop("error! V must be either one of the strings 'alr', 'ilr' or 'clr', or else the matrix Psi of coordinate definition!")
  cn = colnames(COKresult)
  prednames = colnames(COKresult)[grep("pred", cn)]
  varnames = colnames(COKresult)[grep("var", cn)]
  covnames = colnames(COKresult)[grep("cov", cn)]
  coordnames = colnames(COKresult)[-c(grep("cov", cn), grep("var", cn), grep("pred", cn))]
  D = length(prednames)+1
  if(is.null(orignames) & is.matrix(V)) orignames = rownames(V)
  if(is.null(orignames)) orignames = paste("v", 1:D, sep="")
  if(length(orignames)!=D) stop("names provided not consistent with number of logratio variables. Did you forget the rest?")
  prefix = "ilr"
  prediccions = COKresult[,prednames, drop=FALSE]
  if(is.character(V)){
    if(V=="ilr"){
      V = ilrBase(prediccions)
    }else if(V=="alr"){
      V = rbind(diag(D-1), -1)
      prefix = "alr"
    }else if(V=="clr"){
      V = (diag(ncol(prediccions))-matrix(1/D, ncol=D, nrow=D))[, -D]
      prefix = "clr"
    }
  }
  Vsvd = svd(V)
  W = with(Vsvd, v %*% diag(ifelse(d>tol, 1/d, 0)) %*% t(u) )
  colnames(W) = orignames
  if(nscore){
    ## space to back-transform the predictions
    #if(is.null(gg))stop("To apply a nscore backtransformation, the gstat object must be provided!")
    #for(i in 1:(D-1)){
    #  nsc = list()
    #  nsc$trn.table = attr(gg$data[[i]]$data@data[,i],"trn.table")
    #  prediccions[,i] = backtr(scores=prediccions[,i], nscore=nsc,...)
    #}  
    stop("gsi.gstatCokriging2compo: use of 'nscore' is deprecated")
  }
  rg = clrInv(as.matrix(prediccions) %*% W)
  
  # add geographic coordinates as an attribute
  attr(rg,"coords") = COKresult[,coordnames]
  
  if(!nscore){
    # add cokriging variance matrices as an attribute as well
    noms = sub(".pred", "", prednames)
    cvmat = array(0, dim=c(nrow(rg), D-1, D-1), dimnames=list(NULL, noms, noms))
    vrs = COKresult[,varnames, drop=FALSE]
    colnames(vrs) = sub(".var", "", varnames)
    for(ivr in noms){
      cvmat[ ,ivr, ivr] = vrs[,ivr]
    }
    cvs = COKresult[,covnames, drop=FALSE]
    colnames(cvs) = sub("cov.", "", covnames)
    for(ivr in noms){
      for(jvr in noms){
        if(ivr!=jvr){
          dosnoms = c(paste(ivr, jvr, sep="."), paste(jvr, ivr, sep="."))
          quin = dosnoms[dosnoms %in% colnames(cvs)]
          cvmat[ ,ivr, jvr] = cvs[,quin]
        } 
      }
    }
    attr(rg,"krigVar") = cvmat
  }
  class(rg) = c("spatialGridAcomp","acomp")
  return(rg)
}



#' @export
gsi.gstatCokriging2rmult <- function(COKresult, ...) UseMethod("gsi.gstatCokriging2rmult", COKresult)

#' @describeIn gsi.gstatCokriging2compo Reorganisation of cokriged multivariate data
#' @method gsi.gstatCokriging2rmult default
#' @export
gsi.gstatCokriging2rmult.default <- function(COKresult, ...){
  stopifnot(is(COKresult, "Spatial"), "data" %in% slotNames(COKresult))
  coord = sp::coordinates(COKresult)
  res = gsi.gstatCokriging2rmult(cbind(coord, COKresult@data), ...)
  if("grid" %in% slotNames(COKresult)){
    G = ncol(coord)
    dr = 1+(getGridTopology(COKresult)@cellsize<0)
    rp = paste(c("top","bottom")[dr[1]], c("left","right")[dr[2]], sep="")
    if(G==3) rp = paste(rp, c("deep", "surf")[dr[3]], sep="")
    res = setGridOrder(res, refpoint = rp, cycle = 1:G)
  } 
  return(res)
  
}

#' @describeIn gsi.gstatCokriging2compo Reorganisation of cokriged multivariate data
#' @method gsi.gstatCokriging2rmult data.frame
#' @export
gsi.gstatCokriging2rmult.data.frame = function(COKresult, # output of predict.gstat
                                    nscore=FALSE, # were the data NS-transformed?
                                    gg=NULL, # in the case that NS was applied, provide the gstat object!
                                    ... # further arguments needed for nscore
){
  cn = colnames(COKresult)
  prednames = colnames(COKresult)[grep("pred", cn)]
  varnames = colnames(COKresult)[grep("var", cn)]
  covnames = colnames(COKresult)[grep("cov", cn)]
  coordnames = colnames(COKresult)[-c(grep("cov", cn), grep("var", cn), grep("pred", cn))]
  D = length(prednames)
  noms = sub(".pred", "", prednames)
  
  prediccions = rmult(COKresult[,prednames, drop=FALSE])
  if(nscore){
    ## space to back-transform the predictions
    #if(is.null(gg))stop("To apply a nscore backtransformation, the gstat object must be provided!")
    #for(i in 1:(D-1)){
    #  nsc = list()
    #  nsc$trn.table = attr(gg$data[[i]]$data@data[,i],"trn.table")
    #  prediccions[,i] = backtr(scores=prediccions[,i], nscore=nsc,...)
    #}  
    stop("gsi.gstatCokriging2compo: use of 'nscore' is deprecated")
  }
  rg = prediccions
  colnames(rg) = noms  
  # add geographic coordinates as an attribute
  attr(rg,"coords") = COKresult[,coordnames]
  
  if(!nscore){
    # add cokriging variance matrices as an attribute as well
    cvmat = array(0, dim=c(nrow(rg), D, D), dimnames=list(NULL, noms, noms))
    vrs = COKresult[,varnames, drop=FALSE]
    colnames(vrs) = sub(".var", "", varnames)
    for(ivr in noms){
      cvmat[ ,ivr, ivr] = vrs[,ivr]
    }
    cvs = COKresult[,covnames, drop=FALSE]
    colnames(cvs) = sub("cov.", "", covnames)
    for(ivr in noms){
      for(jvr in noms){
        if(ivr!=jvr){
          dosnoms = c(paste(ivr, jvr, sep="."), paste(jvr, ivr, sep="."))
          quin = dosnoms[dosnoms %in% colnames(cvs)]
          cvmat[ ,ivr, jvr] = cvs[,quin]
        } 
      }
    }
    attr(rg,"krigVar") = cvmat
  }
  class(rg) = c("spatialGridRmult","rmult")
  return(rg)
}


#### spatial grids for package "compositions" ------------
# (first attempt: ideally they should be extensions of Spatial)
#' Construct a regionalized composition / reorder compositional simulations
#' 
#' Connect some coordinates to a composition (of hard data, of predictions
#' or of simulations); currently, the coordinates
#' are stored in an attribute and the dataset is given a complex S3 class.
#' This functionality **will** change in the future, to make use of package
#' "sp" classes.  
#'
#' @param coords coordinates of the locations
#' @param compo (observed or predicted) compositional data set; or else array of 
#' simulated compositions
#' @param dimcomp which of the dimensions of `compo` does correspond to the 
#' parts of the compositon? 
#' @param dimsim if `compo` contains simulations, which of its dimensions does 
#' run across the realisations? leave it as NA if `compo` has observations or predictions.
#' @return A (potentially transposed/aperm-ed) matrix of class  c("spatialGridAcomp","acomp")
#' with the coordinates in an extra attribute "coords".
#' @seealso [image_cokriged.spatialGridAcomp()] for an example; [gsi.gstatCokriging2compo()] to
#' restructure the output from [gstat::predict.gstat()] confortably
#' @export
spatialGridAcomp = function(coords, compo, dimcomp=2, dimsim=NA){
  res = compo
  if(is.na(dimsim)){
    if(dimcomp==1){
      res = t(res)
    }
  }else{
    res = aperm(res, order(c(1,dimsim, dimcomp)))
  }
  attr(res,"coords") = coords
  class(res) = c("spatialGridAcomp","acomp")
  return(res)
}


#' Construct a regionalized multivariate data
#' 
#' Connect some coordinates to a multivariate data set (of hard data, of predictions
#' or of simulations); currently, the coordinates
#' are stored in an attribute and the dataset is given a complex S3 class.
#' This functionality **will** change in the future, to make use of package
#' "sp" classes.  
#'
#' @param coords coordinates of the locations
#' @param data (observed or predicted) rmult or matrix data set; or else array of 
#' simulated rmult /real-valued multivariate data
#' @param dimcomp which of the dimensions of `data` does correspond to the 
#' variables? 
#' @param dimsim if `data` contains simulations, which of its dimensions does 
#' run across the realisations? leave it as NA if `data` has observations or predictions.
#' @return A (potentially transposed/aperm-ed) matrix of class  c("spatialGridAcomp","acomp")
#' with the coordinates in an extra attribute "coords".
#' @seealso [image_cokriged.spatialGridRmult()] for an example; [gsi.gstatCokriging2rmult()] to
#' restructure the output from [gstat::predict.gstat()] confortably
#' @export
spatialGridRmult = function(coords, data, dimcomp=2, dimsim=NA){
  res = data
  if(is.na(dimsim)){
    if(dimcomp==1){
      res = t(res)
    }
  }else{
    res = aperm(res, order(c(1,dimsim, dimcomp)))
  }
  attr(res,"coords") = coords
  class(res) = c("spatialGridRmult","rmult")
  return(res)
}


#### image functions for cokriging results ---------
#' Plot an image of gridded data
#' 
#' Plot an image of one variable (possibly, one realisation) of output
#' of cokriging or cosimulation functions.
#'
#' @param x object with the interpolated / simulated data; currently there are methods
#' for "spatialGridAcomp" and "spatialGridRmult", but the default method is able to 
#' deal with "SpatialPointsDataFrame", "SpatialPixelsDataFrame" and "SpatialGridDataFrame"
#' objects, and with the "data.frame" output of [gstat::predict.gstat()] and 
#' [predict_gmSpatialModel]
#' @param ... generic functionality, currently ignored
#'
#' @return Invisibly, a list with elements `breaks` and `col` containing the breaks 
#' and hexadecimal colors finally used for the z-values of the image. Particularly 
#' useful for plotting other plotting elements on the same color scale.
#' @export
#' @importFrom stats quantile 
#' @importFrom graphics layout image title axis abline matplot matlines matpoints
#' 
#' @examples
#' \dontrun{
#' getTellus(cleanup=TRUE, TI=TRUE)
#' load("Tellus_TI.RData")
#' head(Tellus_TI)
#' coords = as.matrix(Tellus_TI[,1:2])
#' compo = compositions::acomp(Tellus_TI[,3:7])
#' dt = spatialGridAcomp(coords=coords, compo=compo)
#' image_cokriged(dt, ivar="MgO") # equi-spaced
#' image_cokriged(dt, ivar="MgO", breaks = NULL) # equi-probable
#' }
image_cokriged <- function(x, ...)  UseMethod("image_cokriged", x)


#' @describeIn image_cokriged Plot an image of gridded data
#' @param ivar which variable do you want to plot?
#' @param breaks either the approximate number of breaks, or the vector of
#'  exact breaks to use for the plotting regions of the chosen variable
#' @param col vector of colors to use for the image
#' @param legendPropSpace which proportion of surface of the device should be used
#' for the legend? trial and error might be necessary to adjust this to your needs
#' @param legendPos where do you want your legend? one of c("top","left","right","bottom")
#' @param main main title for the plot
#' @method image_cokriged default 
#' @export
image_cokriged.default <- function(x, ivar=3, 
                                   breaks=quantile(as.data.frame(x)[,ivar], probs=c(0:10)/10, na.rm=TRUE),
                                   col = spectralcolors(length(breaks)-1),
                                   legendPropSpace = 0.2, legendPos="top",
                                   main=ifelse(is.character(ivar),ivar,colnames(x)[ivar]),
                                   ...){
  if(is(x,"SpatialGridDataFrame")){
    coords = sp::coordinates(x)
    # choose variables
    z = x@data[,ivar]
  }else if(is(x,"SpatialPixelsDataFrame")){
    coords0 = as(x, "SpatialPoints")
    gt = getGridTopology(x)
    idx = getGridIndex(cc=sp::coordinates(coords0), grid=gt)
    # choose variables
    z = rep(NA, prod(gt@cells.dim))
    z[idx] = x@data[,ivar]
    coords = sp::coordinates(gt)
  }else if(is(x,"SpatialPointsDataFrame")){
    coords0 = as(x, "SpatialPoints")
    gt = points2grid(coords0)
    idx = getGridIndex(cc=sp::coordinates(coords0), grid=gt)
    # choose variables
    z = rep(NA, prod(gt@cells.dim))
    z[idx] = x@data[,ivar]
    coords = sp::coordinates(gt)
  }else{
    coords = x[,1:2]
    # choose variables
    z = x[,ivar]
  }
  
  
  if(length(breaks)==1) breaks = pretty(z, n = breaks)
  
  if(length(breaks)-length(col)!=1) col = colorRampPalette(col)(length(breaks)-1)
  
  if(is.logical(legendPos)){
    if(legendPos) legendPos = "top"
  }
  if(!is.logical(legendPos)){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    
    # make sure where the color legend goes
    par(oma=c(1,1,1,1))
    if(legendPos %in% c("top","left")){
      ord = 1:2
      space = c(legendPropSpace,1-legendPropSpace)
    }else if(legendPos %in% c("right","bottom")){
      ord = 2:1
      space = c(1-legendPropSpace,legendPropSpace)
    }else stop("legend must be one of 'top', 'left', 'bottom' or 'right'")
    if(legendPos  %in% c("top","bottom")){
      ord = matrix(ord,ncol=1)
      layout(ord,widths=1, heights=space)
      a = ifelse(legendPos=="top", 4,1)
      par(mar=c(4.5-a,5,a,0.1))
      image(breaks, c(0,1), cbind(breaks, breaks), breaks=breaks, col=col, xlab="", ylab="", 
            yaxt="n", xaxs="i", xlim=range(breaks), oldstyle=TRUE)
      par(mar=c(4.5,5,1+a,0.1))
    }else{
      ord = matrix(ord,nrow=1)
      layout(ord,widths=space, heights=1)
      par(mar=c(4.5,5,2,0.1))
      image(c(0,1), breaks, rbind(breaks, breaks), breaks=breaks, col=col, xlab="", ylab="", xaxt="n",
            yaxs="i", ylim=range(breaks), oldstyle=TRUE)
      par(mar=c(4.5,5,2,0.1))
    }
  }
  # sort the data so that it can be plotted with image 
  oo = order(coords[,2], coords[,1])
  coords = coords[oo,]
  z = z[oo]
  # dimensionalizse 
  xx = sort(unique(coords[,1]))
  yy = sort(unique(coords[,2]))
  dim(z) = c(length(xx),length(yy))
  # plot!
  image(xx,yy,z, breaks=breaks, col=col,
        xlab=colnames(coords)[1],ylab=colnames(coords)[2], asp=1)
  title(main=main, outer=TRUE, line=-1)
  invisible(list(breaks=breaks, col=col))
}



#' @describeIn image_cokriged method for spatialGridRmult objects
#' @param isim in case of simulated output, which simulation?
#' @param mask optional mask object if `x` is of class "spatialGridAcomp" or
#' "spatialGridRmult", and they have been masked  (see [setMask()])
#' @method image_cokriged spatialGridRmult
#' @export
image_cokriged.spatialGridRmult<- function(x, ivar=1, isim=NULL, breaks=10, mask=attr(x, "mask"),
                                             col = spectralcolors(length(breaks)-1),
                                             legendPropSpace = 0.2, legendPos="top",
                                             main=ifelse(is.character(ivar),ivar,dimnames(x)[[length(dimnames(x))]][ivar]),
                                             ...){
    # define data, and if necessary, breaks:
    if(is.null(isim)){
      X = x[,ivar]
    }else{
      X = x[,isim, ivar]
    }
    if(!is.null(mask) & inherits(mask,"mask")){
      if(!is.null(attr(mask, "fullgrid"))){
        X = unmask(data.frame(X), mask=mask)[,1]
        coords = sp::coordinates(attr(mask, "fullgrid"))
      }else stop("mask should be a 'mask' object with stored fullgrid!")
    }else{
      coords = attr(x, "coords")
    }
    oo = order(coords[,2], coords[,1])
    X = X[oo]
    xx = sort(unique(coords[,1]))
    yy = sort(unique(coords[,2]))
    dim(X) = c(length(xx),length(yy))
    
    
    if(is.null(breaks))
      breaks = quantile(X, probs=c(0:10)/10, na.rm=TRUE)

    if(length(breaks)==1) breaks = pretty(X, n = breaks)
    
    if(length(breaks)-length(col)!=1) col = colorRampPalette(col)(length(breaks)-1)
    # make sure where the color legend goes
    par(oma=c(1,1,1,1))
    if(legendPos %in% c("top","left")){
      ord = 1:2
      space = c(legendPropSpace,1-legendPropSpace)
    }else if(legendPos %in% c("right","bottom")){
      ord = 2:1
      space = c(1-legendPropSpace,legendPropSpace)
    }else stop("legend must be one of 'top', 'left', 'bottom' or 'right'")
    if(legendPos  %in% c("top","bottom")){
      ord = matrix(ord,ncol=1)
      layout(ord,widths=1, heights=space)
      a = ifelse(legendPos=="top", 4,1)
      par(mar=c(4.5-a,5,a,0.1))
      image(breaks, c(0,1), cbind(breaks, breaks), breaks=breaks, col=col, xlab="", ylab="", 
            yaxt="n", xaxs="i", xlim=range(breaks))
      par(mar=c(4.5,5,1+a,0.1))
    }else{
      ord = matrix(ord,nrow=1)
      layout(ord,widths=space, heights=1)
      par(mar=c(4.5,5,2,0.1))
      image(c(0,1), breaks, rbind(breaks, breaks), breaks=breaks, col=col, xlab="", ylab="", xaxt="n",
            yaxs="i", ylim=range(breaks))
      par(mar=c(4.5,5,2,0.1))
    }
    image(xx, yy, X, breaks=breaks, col=col, 
          xlab=colnames(coords)[1],
          ylab=colnames(coords)[2], asp=1)
    title(main=main, outer=TRUE, line=-1)
    invisible(list(breaks=breaks, col=col))
  }



#' @describeIn image_cokriged  method for spatialGridAcomp objects
#' @param isim in case of simulated output, which simulation?
#' @param mask optional mask object if `x` is of class "spatialGridAcomp" or
#' "spatialGridRmult", and they have been masked  (see [setMask()])
#' @method image_cokriged spatialGridAcomp
#' @export
image_cokriged.spatialGridAcomp <- image_cokriged.spatialGridRmult


### grid sorting --------


#' Reorder data in a grid
#' 
#' Reorder the data in a compact grid, changing between ordering specifications
#'
#' @param x gridded data
#' @param grid grid topology underlying
#' @param orderIn current ordering description (see [setGridOrder()])
#' @param orderOut desired output ordering description (see [setGridOrder()])
#'
#' @return the data from `x` (typically a matrix), but reordered as `orderOut`
#' @export
#' @seealso [setGridOrder()] for ways of specifying the grid ordering
#' @importClassesFrom sp Spatial
#' @importFrom sp getGridTopology proj4string point.in.polygon 
#' @importFrom sp getGridIndex points2grid
#' @examples
#' \dontrun{
#' getTellus(cleanup=TRUE, TI=TRUE)
#' load("Tellus_TI.RData")
#' coords = as.matrix(Tellus_TI[,1:2])
#' compo = compositions::acomp(Tellus_TI[,3:7])
#' dt = spatialGridAcomp(coords=coords, compo=compo)
#' image_cokriged(dt, ivar="MgO", breaks = NULL) 
#' x = sort(unique(coords[,1]))
#' y = sort(unique(coords[,2]))
#' x0 = c(min(x), min(y))
#' Ax = c(mean(diff(x)), mean(diff(y)))
#' n = c(length(x), length(y))
#' gr = sp::GridTopology(cellcentre.offset=x0, cellsize=Ax, cells.dim=n)
#' dt0 = sortDataInGrid(Tellus_TI, grid=gr, orderIn=gridOrder_array(2), 
#'                     orderOut=list(refpoint="bottomright", cycle=2:1))
#' coords = as.matrix(dt0[,1:2])
#' compo = compositions::acomp(dt0[,3:7])
#' spatialGridAcomp(coords=coords, compo=compo)
#' }
sortDataInGrid = function(x, grid=attr(x, "grid"), 
                          orderIn=attr(x, "gridOrder"), 
                          orderOut=list(refpoint="bottomleft", cycle=1:2)){
  if(is.null(grid)) stop("sortDataInGrid: grid argument cannot be null")
  if(is.null(orderIn)) orderIn = list(refpoint = "topleft", cycle = 1:min(dim(as.data.frame(grid))))
    if(is(grid,"SpatialGrid")) return(sortDataInGrid(x=x, grid=getGridTopology(grid), orderIn = orderIn, orderOut = orderOut)) 
  if(is(x, "Spatial")) return(sortDataInGrid(x=x@data, grid=grid, orderIn = orderIn, orderOut = orderOut))
  if(length(dim(x))==0) x = matrix(x, nrow=length(x), ncol=1)
  # indexes array
  G = length(orderIn$cycle)
  ocycle = (1:G)[orderIn$cycle]
  i = aperm( array(1:nrow(x), dim = grid@cells.dim[orderIn$cycle] ), perm = ocycle)
  # flip 
  dd = dim(i)
  irow = 1:(dd[1])
  icol = 1:(dd[2])
  if(length(grep("right",orderIn$refpoint))>0) irow = rev(irow)
  if(length(grep("top",orderIn$refpoint))>0) icol = rev(icol)
  if(grid@cellsize[1]<0) irow = rev(irow)
  if(grid@cellsize[2]<0) icol = rev(icol)
  if(G==3){
    ideep = 1:(dd[3])
    if(length(grep("surf",orderIn$refpoint))>0) ideep = rev(ideep)
    if(grid@cellsize[3]<0) ideep = rev(ideep)
    i = i[irow, icol, ideep]
  }else{
    i = i[irow, icol]
  }
  y = x[i,, drop=T]
  if(all(orderOut$cycle==1:G) & orderOut$refpoint=="bottomleft") return(y)

  ocycle = (1:G)[orderOut$cycle]
  j = aperm( array((1:length(i))[i], dim = grid@cells.dim[orderOut$cycle] ), perm = ocycle)
  
  dd = dim(j)
  jrow = 1:(dd[1])
  jcol = 1:(dd[2])
  if(length(grep("right",orderOut$refpoint))>0) jrow = rev(jrow)
  if(length(grep("top",orderOut$refpoint))>0) jcol = rev(jcol)
  if(G==3){
    jdeep = 1:(dd[3])
    if(length(grep("surf",orderOut$refpoint))>0) jdeep = rev(jdeep)
    j = j[jrow, jcol, jdeep]
  }else{
    j = j[jrow, jcol]
  }
  y = x[j,, drop=T]
  return(y)
}


