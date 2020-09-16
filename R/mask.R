### Constructors for masks --------------
#' Constructs a mask for a grid
#'
#' @param grid a grid, see details for more info
#' @param method which construction method? currently one of 'maxdist', 'sillprop' or 'point2polygon'
#' @param maxval for maxdist and sillprop methods, maximum reference value
#' @param x extra information for the grid construction, see details
#'
#' @return a logical vector with as many elements as points in the grid, with TRUE for
#' those points within the mask, and FALSE for those outside the mask.
#' @details Method 'maxdist' defines the mask as all points within a maximum distance 
#' (must be given in \code{maxval}) from the reference data (given in \code{x}: this is expected 
#' to be the original complete data, with coordinates and variables). For method 'sillprop'
#' the mask is defined by those points which total kriging variance is below
#' a fixed proportion (given in \code{maxval}, default=0.99) of the total variogram
#' model sill (variogram model given in \code{x}, of class "variogramModelList"). 
#' In this method, the argument \code{grid} is expected to be the output of a cokriging
#' analysis. Finally, method 'point2poly' created the mask by taking the points internal
#' to a "SpatialPolygon" object (given in \code{x}).   
#' @export
#' @family masking functions
#'
#' @examples
#' ## with data.frame
#' x = 1:23
#' y = 1:29
#' xy = expand.grid(x=x, y=y)
#' xyz.df = data.frame(xy, z = rnorm(29*23)*ifelse(abs(xy$x-xy$y)<3, 1, NA)+(xy$x+xy$y)/2)
#' mask.df = constructMask(grid = xy, method = "maxdist", maxval = 3, x=xyz.df)
#' image(mask.df)
#' par(mfrow=c(1,1))
#' mask.df
#' xyz.df.masked = setMask(xyz.df, mask.df)
#' dim(xyz.df.masked)
#' summary(xyz.df.masked)
#' xyz.df.unmasked = unmask(xyz.df.masked)
#' dim(xyz.df.unmasked)
#' length(x)*length(y)
#' summary(xyz.df.unmasked)
#' ## with SpatialGrid
#' library(sp)
#' library(magrittr)
#' xy.sp = sp::SpatialPoints(coords = xy)
#' meandiff = function(x) mean(diff(x))
#' xy.gt = GridTopology(c(min(x),min(y)), c(meandiff(x), meandiff(y)), c(length(x),length(y)))
#' aux = sp::SpatialPixelsDataFrame(grid = xy.gt, data=xyz.df, points = xy.sp)
#' xyz.sgdf = as(aux, "SpatialGridDataFrame")
#' image_cokriged(xyz.sgdf, ivar="z")
#' par(mfrow=c(1,1))
#' ms = function(x) sortDataInGrid(x, grid=xy.gt)
#' mask.gt = constructMask(grid = xy.gt, method = "maxdist", maxval = 3, x=xyz.sgdf)
#' image(x,y,matrix(ms(xyz.sgdf@data$z), nrow=23, ncol=29))  
#' image(x,y,matrix(ms(mask.gt), nrow=23, ncol=29))  
#' image(mask.gt)
#' par(mfrow=c(1,1))
#' xyz.sgdf.masked = setMask(x = xyz.sgdf, mask = mask.gt)
#' getMask(xyz.sgdf.masked)
#' image(x,y,matrix(ms(xyz.sgdf@data$z), nrow=23, ncol=29))  
#' points(xyz.sgdf.masked@coords)
constructMask = function(grid, method="maxdist", maxval=NULL, x=NULL){
  methods = c("maxdist", "sillprop", "point2poly") 
  m = methods[pmatch(method, methods)] 
  if(is(grid,"GridTopology")){
    grid0 = sp::coordinates(sp::SpatialGrid(grid))
  }else if(is(grid,"Spatial") & ("data" %in% slotNames(grid))){
    grid0 = data.frame(sp::coordinates(grid), grid@data)
    if(is.null(x)) x = grid0
  }else if(is(grid,"Spatial")){
    grid0 = sp::coordinates(grid)
  }else if(is.data.frame(grid)){
    grid0 = grid
    if(is.null(x)) x = grid0
  }else stop("constructMask: object 'grid' could not be interpreted")
  if(is.na(m))
    stop('constructMask: method should be one of c("maxdist", "sillprop", "point2poly")')
  if(m=="maxdist"){
    if(is.null(maxval))
      stop("constructMask: maxdist method requires a maxval=maximum distance to location")
    if(is(x, "SpatialPointsDataFrame")){
      x = data.frame(sp::coordinates(x), x@data)
    }
    x = try(as.data.frame(x))
    if(class(x)=="try-error") stop("constructMask: provided object x should be a data.frame or convertible to it for method 'maxdist'")
    out = gsi.masking.nearest(grid0, x, maxdist=maxval)
  }else if(m=="sillprop"){
    if(is.null(x))
      stop("constructMask: sillprop method requires a variogram model")
    if(class(x)=="gstat") x = x$model
    if(is(x,"gmSpatialModel")) x = x@model@structure
    if(is(x, "ModelStructuralFunctionSpecification")) as.variogramModel(x)
    maxval = ifelse(is.null(maxval), 0.99, maxval)
    out = gsi.masking.cokriged(out.ck=grid0, vgmodel=x, sillFraction=maxval)
  }else if(m=="point2polygon"){
    out = gsi.masking.polygon(grid0, x)
  }else
    stop("constructMask: method must be one of 'maxdist' or 'sillprop'")
  out[is.na(out)]=FALSE
  attr(out, "fullgrid") = grid
  class(out) = "mask"
  return(out)
}


#' Image method for mask objects
#' 
#' Plot a mask of a 2D grid; see [constructMask()] for an example
#'
#' @param x a mask
#' @param col a two-color vector for the color (oustide, inside) the mask
#' @param ... ignored
#'
#' @return nothing
#' @export
image.mask <- function(x, col=c(NA,2), ...){
  grid = attr(x, "fullgrid")
  if(is(grid,"GridTopology")){
    grid0 = sp::coordinates(sp::SpatialGrid(grid))
    o = order(-grid0[,2],+grid0[,1])
    grid0 = grid0[o,]
  }else if(is(grid,"SpatialGrid")){
    grid0 = sp::coordinates(grid)
    o = order(-grid0[,2],+grid0[,1])
    grid0 = grid0[o,]
  }else if(is(grid,"Spatial")){
    grid0 = sp::coordinates(grid)
  }else if(is.data.frame(grid)){
    grid0 = grid
  }
  image_cokriged(cbind(grid0, mask=unclass(x)), ivar="mask", breaks = c(-0.0001, mean(unclass(x)), 1.0001),
                 col = col)
}

### masks all predictions which total variance larger than a certain trace sill fraction 
gsi.masking.cokriged = function(out.ck, vgmodel, sillFraction=0.99){
  idpreds = grep( "pred", colnames(out.ck))
  idvars = grep("var", colnames(out.ck))
  variables = sub(".var","",colnames(out.ck)[idvars])
  if(length(idvars)==0) stop("method 'sillprop' requires the output of a (co)kriging")
  maxvar = sum(sapply(variables, function(i) sum(vgmodel[[i]]$psill)))
  tracevar = rowSums(out.ck[,idvars])
  mask = tracevar<=(sillFraction*maxvar)
  return(mask)
} 


gsi.masking.nearest = function(grid, x, maxdist){
  if(!requireNamespace("FNN", quietly = TRUE)) stop("constructMask: method='maxdist' requires package 'KNN' installed")
  x = as.matrix(x)
  coordnames = intersect(colnames(grid), colnames(x))
  varnames = setdiff(colnames(x), coordnames)
  aux = matrix(NA, nrow=nrow(grid), ncol = length(varnames))
  colnames(aux) = varnames
  expanded.grid <- cbind(grid, aux)
  aux3 = FNN::get.knnx(expanded.grid[,coordnames], x[,coordnames], k=1, algo="kd_tree")
  expanded.grid[aux3$nn.index, varnames] = x[,varnames]
  expanded.grid = cbind(expanded.grid,0)
  colnames(expanded.grid)[ncol(expanded.grid)] = "Mask"
  expanded.grid[complete.cases(expanded.grid),ncol(expanded.grid)] = 1
  expanded.grid[FNN::get.knnx(expanded.grid[complete.cases(expanded.grid),coordnames],
                         expanded.grid[,coordnames], k=1, algo="kd_tree")$nn.dist<=maxdist,
                ncol(expanded.grid)] = 1
  return( as.logical(expanded.grid[,"Mask"]) )
}


gsi.masking.polygon = function(grid, poly){
  requireNamespace("sp", quietly = TRUE)
  poly = try(as(poly, "SpatialPolygons"))
  if(class(poly)=="try-error")
    stop("object 'poly' cannot be coerced to SpatialPolygons")
  FUN = function(i){
    poly = poly@polygons[[i]]@Polygons[[1]]@coords
    quins = point.in.polygon(grid[,1], grid[,2], poly[,1], poly[,2]) ==1
    return(quins)
  }
  erg = sapply(1:length(poly@polygons), FUN)
  return(apply(erg,1,any))
}


#### getters and setters for masks ---------

#' Get the mask info out of a spatial data object
#' 
#' Retrieve the mask information from an object (if present). See [constructMask()]
#' for examples.
#'
#' @param x a masked object 
#' @return The retrieved mask information from `x`, an object of class "mask" 
#' @export
#' @family masking functions
getMask = function(x) UseMethod(generic = "getMask", x)

#' @describeIn getMask Get the mask info out of a spatial data object
#' @method getMask default
#' @export
getMask.default = function(x) attr(x, "mask")

#' @describeIn getMask Get the mask info out of a SpatialPixelsDataFrame data object
#' @method getMask SpatialPixelsDataFrame
#' @export
#' @importClassesFrom sp SpatialPixelsDataFrame
getMask.SpatialPixelsDataFrame <- function(x){
  coords = sp::coordinates(sp::SpatialGrid(sp::getGridTopology(x)))
  mask = rep(FALSE, nrow(coords))
  mask[x@grid.index] = TRUE
  attr(mask, "fullgrid") = getGridTopology(x)
  class(mask) = "mask" 
  return(mask)
}

#' @describeIn getMask Get the mask info out of a SpatialPixels object
#' @method getMask SpatialPixels
#' @export
#' @importClassesFrom sp SpatialPixels
getMask.SpatialPixels <- getMask.SpatialPixelsDataFrame

#' @describeIn getMask Get the mask info out of a SpatialPointsDataFrame data object
#' @method getMask SpatialPointsDataFrame
#' @export
getMask.SpatialPointsDataFrame = function(x) attr(x@data, "mask")

#' Print method for mask objects
#' 
#' Print method for mask objects. See [constructMask()] for examples.
#' If you want to see the whole content of the mask, then use `unclass(...)`
#'
#' @param x mask to print
#' @param ... ignored
#'
#' @return the summary of number of nodes inside/outside the mask
#' @export
#' @family masking functions
print.mask <- function(x,...){
  print("mask active")
  print(summary(x))
} 



#' Set a mask on an object
#' 
#' Set a mask on an object See [constructMask()] for examples on how to construct masks.
#'
#' @param x an object to mask (for set) or masked (for get)
#' @param ... extra arguments for generic compatibility
#'
#' @return The object `x` appropriately masked (for the setter methods). 
#' @export
#' @family masking functions
setMask <- function(x,...) UseMethod("setMask", x)


#' @describeIn setMask Set a mask on an object
#' @export
#' @method setMask default
#' @param mask the mask to impose on `x`
#' @param coordinates for some of the methods, it is important to specify the names or indices
#' of the columns containing the geographic coordinates (only `setMask.data.frame`) or else
#' to specify the matrix of spatial coordinates (all `setMask` methods including it)
setMask.default <- function(x, mask, coordinates = 1:2, ...){
  x = as.data.frame(x)
  if(class(mask)=="mask") attributes(mask) = NULL
  if(is.null(dim(coordinates) )){
    fullgrid = x[,coordinates]
  }else{
    fullgrid = coordinates
  }
  outdata = x[mask,,drop=FALSE]
  attr(mask, "fullgrid") = fullgrid
  attr(outdata, "mask") = mask
  return(outdata)
}

#' @describeIn setMask Set a mask on a data.frame object
#' @method setMask data.frame
#' @export
setMask.data.frame <- setMask.default

#' @describeIn setMask Set a mask on a DataFrameStack object
#' @method setMask DataFrameStack
#' @export
setMask.DataFrameStack <- function(x, mask, coordinates=attr(x, "coordinates"), ...){
  if(class(mask)=="mask") attributes(mask) = NULL
  cc = coordinates
  x = x[mask,,drop=FALSE]
  attr(mask, "fullgrid") = cc
  attr(x, "mask") = mask
  return(x)
}

#' @describeIn setMask Set a mask on a SpatialGrid object
#' @method setMask SpatialGrid
#' @export
#' @importClassesFrom sp SpatialGrid
setMask.SpatialGrid <- function(x, mask, ...){
  cc = sp::coordinates(x)
  r = order(+cc[,2],+cc[,1])
  o = 1:nrow(cc)
  r = o[r]
  if(class(mask)=="mask") attributes(mask) = NULL
  maskaux = mask[o]
  cc = cc[maskaux,, drop=FALSE]
  cc = sp::SpatialPoints(coords = cc, proj4string = sp::CRS(sp::proj4string(x)), 
                         bbox = sp::bbox(x))
  if("data" %in% slotNames(x)){
    dt = x@data
    dt = dt[maskaux,, drop=FALSE]
    erg = sp::SpatialPixelsDataFrame(points=cc, data = dt, 
                                     grid = sp::getGridTopology(x), 
                                     proj4string =sp::CRS(sp::proj4string(x)))
  }else{
    erg = sp::SpatialPixels(points = cc, 
                            grid = sp::getGridTopology(x), 
                            proj4string = sp::CRS(sp::proj4string(x) ) )
  }
  return(erg)
}

#' @describeIn setMask Set a mask on a GridTopology object
#' @method setMask GridTopology
#' @export
#' @importClassesFrom sp GridTopology
setMask.GridTopology <- function(x, mask, ...){
  setMask(sp::SpatialGrid(x), mask, ...)
}

#' @describeIn setMask Set a mask on a SpatialPoints object
#' @method setMask SpatialPoints
#' @export
#' @importClassesFrom sp SpatialPoints
setMask.SpatialPoints <- function(x, mask, ...){
  cc = sp::coordinates(x)
  if(class(mask)=="mask") attributes(mask) = NULL
  cc = cc[mask,,drop=FALSE]
  if("data" %in% slotNames(x)){
    dt = x@data
    dt = dt[mask,,drop=FALSE]
    attr(mask, "fullgrid") = sp::coordinates(x)
    attr(dt, "mask") = mask
    erg = sp::SpatialPointsDataFrame(coords = cc, data = dt, bbox = sp::bbox(x), 
                                 proj4string = sp::CRS(sp::proj4string(x)) )
  }else{
    erg = sp::SpatialPoints(coords = cc, bbox = sp::bbox(x), 
                        proj4string = sp::CRS(sp::proj4string(x)) )
  }
  return(erg)
}


#### unmasking function ---------

#' @describeIn unmask.data.frame Unmask a masked object
#' @param ... arguments for generic functionality
#' @export
unmask <- function(x,...) UseMethod("unmask", x)

#' Unmask a masked object
#' 
#' Unmask a masked object, i.e. recover the original grid and extend potential
#' data containers associated to it with NAs. See examples in [constructMask()]
#'
#' @param x a masked object
#' @param mask the mask; typically has good defaults
#' @param fullgrid the full grid; typically has good defaults
#' @param forceCheck if `fullgrid` is provided, should the coordinates provided 
#' in `x` and in `fullgrid` be cross-checked to ensure that they are given in 
#' compatible orders? See [sortDataInGrid()] and [setGridOrder()] for controlling
#' the ordering of vectors and grids.
#' @family masking functions
#' @method unmask data.frame
#'
#' @return The original grid data and extend potential
#' data containers associated to it with NAs. See examples in [constructMask()].
#' The nature of the output depends on the nature of `x`:
#' a "data.frame" produced a "data.frame";
#' a "unmask.DataFrameStack" produces a "unmask.DataFrameStack"; 
#' a "SpatialPoints" produces a "SpatialPoints"; and finally
#' a "SpatialPixels" produces either a "SpatialPixels" or a "SpatialGrid" (if it is full).
#' Note that only in the case that `class(x)=="SpatialPixels"` is `mask` required,
#' for the other methods all arguments have reasonable defaults.
#' @export
unmask.data.frame <- function(x, mask=attr(x,"mask"), fullgrid = attr(mask, "fullgrid"), 
                              forceCheck=is(fullgrid, "GridTopology"), ...){
  if(is(fullgrid, "GridTopology")) fullgrid = sp::SpatialGrid(fullgrid)
  if(is(fullgrid, "Spatial")) fullgrid = data.frame(sp::coordinates(fullgrid))
  out = data.frame(matrix(NA, ncol=ncol(x)-ncol(fullgrid), nrow=length(mask)))
  out = cbind(fullgrid, out)
  colnames(out) = colnames(x)
  out[mask, colnames(x)] = x
  return(out)
}

#' @describeIn unmask.data.frame Unmask a masked object
#' @method unmask DataFrameStack 
#' @export
unmask.DataFrameStack <- function(x, mask=attr(x,"mask"), fullgrid = attr(mask, "fullgrid"), 
                                  forceCheck=is(fullgrid, "GridTopology"), ...){
  if(is(fullgrid, "GridTopology")) fullgrid = sp::SpatialGrid(fullgrid)
  if(is(fullgrid, "Spatial")) fullgrid = data.frame(sp::coordinates(fullgrid))
  out = data.frame(matrix(NA, ncol=ncol(x), nrow=length(mask)))
  colnames(out) = colnames(x)
  out[mask,] = x
  odimnames = attr(x, "Dimnames")
  rwn <- rownames(fullgrid)
  if(is.null(rwn)) rwn = 1:nrow(fullgrid)
  odimnames = list(rwn, odimnames[[2]], odimnames[[3]])
  names(odimnames) = names(attr(x, "Dimnames"))
  rownames(out) <- rwn
  attr(out, "Dimnames") = odimnames
  attr(out, "stackDim") = attr(x, "stackDim") 
  attr(out, "fullgrid") = fullgrid
  class(out) = class(x)  
  return(out)
}

#' @describeIn unmask.data.frame Unmask a masked object
#' @method unmask SpatialPixels
#' @export
#' @importClassesFrom sp SpatialPixels
unmask.SpatialPixels <- function(x, mask=NULL, fullgrid =attr(mask, "fullgrid"), 
                                 forceCheck=FALSE, ...){
  # store grid topology of the original data
  gtin = sp::getGridTopology(x)
  # extract/construct a grid topology of the fullgrid
  if(is.null(fullgrid)){
    fullgrid = gtin  #... the same as the topology of x if fullgrid absent
  }else if(is.data.frame(fullgrid)){ # or construct a grid if coordinates are provided as data.frame or SpatialPoints(DataFrame)
    fullgrid = sp::SpatialPixels(points=fullgrid, proj4string = sp::CRS(sp::proj4string(x)) ) 
  }else if(is(fullgrid, "SpatialPoints")){
    fullgrid = sp::SpatialPixels(points=sp::coordinates(fullgrid), proj4string = sp::CRS(sp::proj4string(x)) )
  }
  if(is(fullgrid, "SpatialGrid")) fullgrid = sp::getGridTopology(fullgrid)
  # compute number of points
  npoints <- try( prod(fullgrid@cells.dim))
  if(class(npoints)=="try-error") stop("unmask.SpatialPixels: provided fullgrid could not be intepreted as a grid")
  # construct mask
  if(is.null(mask)){
    mask = rep(FALSE, npoints)
    mask[x@grid.index]=TRUE
  } 
  if("data" %in% slotNames(x)){
    # deal with SpatialPixelsDataFrame
    dt = x@data
    X = as.data.frame(matrix(NA, ncol=ncol(dt), nrow=npoints))
    colnames(X) = colnames(dt)
    X[mask,] = dt
    erg = sp::SpatialGridDataFrame(grid=fullgrid, data=X, 
                                   proj4string = sp::CRS(sp::proj4string(x)))
  }else{
    # deal with SpatialPixels
    erg = sp::SpatialGrid(grid=fullgrid, proj4string = sp::CRS(sp::proj4string(x)))
  }
  return(erg)
}


#' @describeIn unmask.data.frame Unmask a masked object
#' @method unmask SpatialPoints
#' @importClassesFrom sp SpatialPoints
unmask.SpatialPoints <- function(x, mask=attr(x@data,"mask"), 
                                 fullgrid = attr(mask, "fullgrid"), 
                                 forceCheck=FALSE, ...){
  # stop("unmask.SpatialPoints: not yet implemented")
  if(is(fullgrid, "GridTopology")){
    # stop("unmask.SpatialPoints: not yet implemented for fullgrid GridTopology")
    
  }else if(is(fullgrid, "SpatialGrid")){
    # stop("unmask.SpatialPoints: not yet implemented for fullgrid SpatialGrid")
    
  }else if(is(fullgrid, "SpatialPoints")){
    dt.x = cbind(sp::coordinates(x))
    if("data" %in% slotNames(x)) dt.x = cbind(dt.x, x@data)
    dt.f = data.frame(sp::coordinates(fullgrid))
    erg = unmask(dt.x, mask=mask, fullgrid=dt.f, forceCheck = forceCheck | is(fullgrid, "SpatialPixels"))
    res = sp::SpatialPoints(coords =erg[,1:ncol(dt.f)], 
                            bbox = sp::bbox(fullgrid), 
                            proj4string = sp::CRS(sp::proj4string(x)))
    if("data" %in% slotNames(x)) res = sp::SpatialPointsDataFrame(coords = res, data=erg[, -(1:ncol(dt.f)), drop=F])  
    return(res)
  }
  warning("unmask.SpatialPoints: strange fullgrid provided; attempting a patch")
  unmask(as.data.frame(x), mask=mask, fullgrid=fullgrid, forceCheck=forceCheck)
}






# 
# ### tests -----------
# if(!exists("do.test")) do.test=FALSE
# if(do.test){
#   
#   ## case data.frame ---
#   # create setting
#   x = 1:23
#   y = 1:29
#   xy = expand.grid(x=x, y=y)
#   xyz.df = data.frame(xy, z = rnorm(29*23)*ifelse(abs(xy$x-xy$y)<3, 1, NA)+(xy$x+xy$y)/2)
#   # image(x,y,matrix(xyz.df$z, nrow=23, ncol=29))  
#   # mask
#   mask.df = constructMask(grid = xy, method = "maxdist", maxval = 3, x=xyz.df)
#   image(mask.df)
#   par(mfrow=c(1,1))
#   image(x,y,matrix(xyz.df$z, nrow=23, ncol=29))  
#   xyz.df.masked = setMask(x = xyz.df, mask = mask.df)
#   points(xyz.df.masked[,1:2])
#   # "interpolate"
#   xyz.df.masked$z <- with(xyz.df.masked, ifelse(is.na(z), (x+y)/2, z) )
#   # unmask
#   xyz.df.unmasked = unmask(xyz.df.masked, mask=mask.df)
#   image(x,y,matrix(xyz.df.unmasked$z, nrow=23, ncol=29))  
#   
#   ## case SpatialPoints ---
#   # create setting
#   xy.sp = SpatialPoints(coords = xy)
#   xyz.spdf = SpatialPointsDataFrame(coords = xy.sp, data=xyz.df)
#   image_cokriged(xyz.spdf, ivar="z")  
#   par(mfrow=c(1,1))
#   # mask
#   mask.sp = constructMask(grid = xy.sp, method = "maxdist", maxval = 3, x=xyz.spdf)
#   image(x,y,matrix(xyz.spdf@data$z, nrow=23, ncol=29))  
#   image(x,y,matrix(mask.sp, nrow=23, ncol=29))  
#   image(mask.sp)
#   par(mfrow=c(1,1))
#   xyz.sp.masked = setMask(x = xyz.spdf, mask = mask.sp)
#   image(x,y,matrix(xyz.spdf@data$z, nrow=23, ncol=29))  
#   points(xyz.sp.masked@coords)
#   # "interpolate"
#   xyz.sp.masked@data$z <- with(xyz.sp.masked@data, ifelse(is.na(z), (x+y)/2, z) )
#   xyz.sp.masked@data <- data.frame(z=xyz.sp.masked@data$z)
#   # unmask
#   xyz.sp.unmasked = unmask(xyz.sp.masked, mask=mask.sp)
#   image(x,y,matrix(xyz.sp.unmasked@data$z, nrow=23, ncol=29)) 
#   image_cokriged(xyz.sp.unmasked, ivar="z")
#   par(mfrow=c(1,1))
#   
#   ## case SpatialGrid ---
#   # create setting
#   meandiff = function(x) mean(diff(x))
#   xy.gt = GridTopology(c(min(x),min(y)), c(meandiff(x), meandiff(y)), c(length(x),length(y)))
#   xyz.sgdf = SpatialPixelsDataFrame(grid = xy.gt, data=xyz.df, points = xy.sp) %>% as("SpatialGridDataFrame")
#   image_cokriged(xyz.spdf, ivar="z")  
#   par(mfrow=c(1,1))
#   # mask
#   ms = function(x) sortDataInGrid(x, grid=xy.gt)
#   mask.gt = constructMask(grid = xy.gt, method = "maxdist", maxval = 3, x=xyz.sgdf)
#   image(x,y,matrix(ms(xyz.sgdf@data$z), nrow=23, ncol=29))  
#   image(x,y,matrix(ms(mask.gt), nrow=23, ncol=29))  
#   image(mask.gt)
#   par(mfrow=c(1,1))
#   xyz.sgdf.masked = setMask(x = xyz.sgdf, mask = mask.gt)
#   image(x,y,matrix(ms(xyz.sgdf@data$z), nrow=23, ncol=29))  
#   points(xyz.sgdf.masked@coords)
#   # "interpolate"
#   z <- with(xyz.sgdf.masked@data, ifelse(is.na(z), (x+y)/2, z) )
#   xyz.sgdf.masked = SpatialPixelsDataFrame(points = xyz.sgdf.masked@coords, 
#                                            data = data.frame(z=z),grid =xy.gt )
#   # image(xyz.sgdf.masked)  
#   # unmask
#   xyz.sp.unmasked = unmask(xyz.sgdf.masked, mask=mask.gt)
#   image(x,y,as.matrix(xyz.sgdf.masked)[,(29:1)])  # logic, but useless
#   image_cokriged(xyz.sp.unmasked, ivar="z")
#   dev.off()
#   
#   ## case GridTopology ---
#   # check 
#   xyz.sg.masked = setMask(x = xy.gt, mask = mask.gt)
#   par(mfrow=c(1,1))
#   plot(xyz.sg.masked)
#   plot(sp::coordinates(xyz.sg.masked))
#   points(sp::coordinates(xyz.sp.masked), pch=4, col=2)
#   
#   ## case DataFrameStack ---
#   xyz.dfs = lapply(1:5, function(i) data.frame(z = rnorm(29*23)*ifelse(abs(xy$x-xy$y)<3, 1, NA)+(xy$x+xy$y)/2) )
#   names(xyz.dfs) = LETTERS[1:5]
#   xyz.dfs = DataFrameStack(xyz.dfs, stackDimName = "sim")
#   # image(x,y,matrix(xyz.df$z, nrow=23, ncol=29))  
#   # mask
#   mask.dfs = constructMask(grid = xy, method = "maxdist", maxval = 3, x=cbind(xy,getStackElement(xyz.dfs,1)))
#   image(mask.dfs)
#   par(mfrow=c(1,1))
#   image(x,y, matrix(unlist(getStackElement(xyz.dfs,3)), nrow=23, ncol=29))  
#   xyz.dfs.masked = setMask(x = xyz.dfs, mask = mask.dfs, coordinates = xy)
#   xy.dfs.masked = xyz.dfs.masked %>% attr("mask") %>% attr("fullgrid") %>% setMask(mask = mask.dfs)
#   xy.dfs.masked %>% points
#   # "interpolate"
#   for(i in dimnames(xyz.dfs.masked)[[stackDim(xyz.dfs.masked)]]){
#     newz = with(data.frame(xy.dfs.masked, getStackElement(xyz.dfs.masked,i)), ifelse(is.na(z), (x+y)/2, z) )
#     xyz.dfs.masked %<>%  setStackElement(i, newz) 
#   }
#   # unmask
#   xyz.dfs.unmasked = unmask(xyz.dfs.masked, mask=mask.dfs)
#   image(x,y,matrix(as.matrix(getStackElement(xyz.dfs.unmasked,4)), nrow=23, ncol=29))  
# }



