# Class DataFrameStack
gsi.dataframestack <- function(x, stackDim=2, dim=NULL, Dimnames=NULL,...){
  l = list(...)
  if(is.null(Dimnames) & "dimnames" %in% names(l)) Dimnames= l$dimnames
  if(is.null(dim)&is.null(Dimnames))
    stop("DataFrameStack: one of dim or Dimnames must be provided")
  if(!is.null(dim)&!is.null(Dimnames)){
    warning("DataFrameStack: dim and dimnames were provided, dim is ignored")
    dim2 = sapply(Dimnames,length)
    dim[dim2!=0] = dim2[dim2!=0]
  }
  dn = sapply(Dimnames, is.null)
  if(any(dn)){
    if(is.null(dim)) stop("DataFrameStack: if Dimnames are not complete, dim is necessary")
    dimnames2 = lapply(dim, function(i) 1:i)
    Dimnames[dn] = dimnames2[dn]
  }else if(length(dn)==0){
    if(is.null(dim)) stop("DataFrameStack: if Dimnames are not provided, dim is necessary")
    Dimnames = lapply(dim, function(i) 1:i)
  }
  if(is.null(dim))
    dim = sapply(Dimnames, length)
  if(!gsi.is.stackDim.defined(stackDim, names(Dimnames), length(dim)))
    stop("DataFrameStack: stacking dimension unidentified")
  attr(x, "Dimnames") = Dimnames
  attr(x, "stackDim") = stackDim  
  class(x) = unique(c("DataFrameStack", class(x)))
  return(x)
}

gsi.is.stackDim.defined = function(stackdim, n=NULL, l=NULL){
  if(is.character(stackdim)){
    return(stackdim %in% n)
  }else{
    return(stackdim>0 & stackdim<(l+1))
  }
}


#' @describeIn DataFrameStack.data.frame create a DataFrameStack from an array
#' @export
DataFrameStack <- function(x,...) UseMethod("DataFrameStack",x)


#' @describeIn DataFrameStack.data.frame create a DataFrameStack from an array
#' @export
as.DataFrameStack <- DataFrameStack

#' Create a data frame stack
#' 
#' Make a stacked data frame, that is, a stack of data.frames representing e.g. repeated measurements,
#' parallel time series, or a stack of multivariate realisation of some random process. If a data frame
#' is analogous to a matrix, a DataFrameStack is analogous to an array. It is highly recommendable to work
#' with named dimensions in stacked data frames.
#'
#' @param x object containing the individual data sets, this can be an arra, a list or a data.frame
#' @param stackDim for "array" and "data.frame" input, which dimension 
#' (name or index) identifies the stacking dimension;
#' this is typically the replication, realisation or time slice index 
#' @param dim for "data.frame" input,  how is the data provided to be arranged in slices?
#' @param Dimnames for "list" and "data.frame input, which names do you want to give 
#' to the resulting array-like structure; note that for input "array" it is necessary 
#' that these names are already given to the input array beforehand!
#' @param ... further parameters, ignored but necessary for generic functionality
#'
#' @return The data provided reorganised as a DataFrameStack, with several additional attributes.
#' @export
#' @method DataFrameStack data.frame
#' @seealso [stackDim()] to get or set the stacking dimension; 
#' [noStackDim()] to get (not set) the non-stacking dimension;
#' [as.array.DataFrameStack()] and [as.list.DataFrameStack()] to
#' convert the stack to an array or a list;
#' [dimnames.DataFrameStack()] to get the dimension names;
#' [`[.DataFrameStack`] to extract rows of a stack;
#' [getStackElement()] and `setStackElement` (same page as `getStackElement`) 
#' to extract/modify an 
#' element of the stack; [gmApply()] to
#' apply any function to the stack, typically element-wise;
#' and [swarmPlot()] to combine plot elements for each stack element
#' into a single plot. 
#' @examples
#' ll = lapply(1:3, function(i) as.data.frame(matrix(1:10+(i-1)*10, ncol=2)))
#' dfs = DataFrameStack(ll, stackDimName="rep")
#' dimnames(dfs)
#' df = as.data.frame(matrix(1:30, ncol=6))
#' dfs = DataFrameStack(df, dimnames = list(obs=1:5, vars=c("A","B"), rep=1:3), stackDim = "rep")
#' dimnames(dfs)
#' ar = array(1:30, dim = c(5,2,3), dimnames=list(obs=1:5, vars=c("A","B"), rep=1:3))
#' dfs = DataFrameStack(ar, stackDim="rep")
#' dimnames(dfs)
DataFrameStack.data.frame <-  gsi.dataframestack

#' @describeIn DataFrameStack.data.frame create a DataFrameStack from an array
#' @method as.DataFrameStack data.frame
#' @export
as.DataFrameStack.data.frame <- gsi.dataframestack


#' @describeIn DataFrameStack.data.frame create a DataFrameStack from an array
#' @method DataFrameStack array
#' @export
DataFrameStack.array <- function(x, stackDim=2, ...){
  if(!gsi.is.stackDim.defined(stackDim, names(dimnames(x)), length(dim(x))))
    stop("DataFrameStack: stacking dimension unidentified")
  aux = x
  dim(x) = c(dim(x)[1], prod(dim(x)[-1]))
  x = data.frame(x)
  attr(x, "Dimnames") = dimnames(aux)
  attr(x, "stackDim") = stackDim 
  class(x) = unique(c("DataFrameStack", "data.frame"))
  return(x)
}

#' @describeIn DataFrameStack.data.frame create a DataFrameStack from an array
#' @method as.DataFrameStack array
#' @export
as.DataFrameStack.array <- DataFrameStack.array


#' @describeIn DataFrameStack.data.frame create a DataFrameStack from a list
#' @param stackDimName for "list" input, which name or index do you want to 
#' give to the stacking dimension?  
#' @export
#' @method DataFrameStack list 
DataFrameStack.list <- function(x, stackDimName = NULL, Dimnames=NULL, ...){
  aux = x[[1]]
  if(length(dim(aux))!=2)
    dim(aux) = c(length(aux),1)
  if(is.null(Dimnames)){
    nm = if(is.null(names(x))) 1:length(x) else names(x)
    rn = if(is.null(rownames(aux))) 1:nrow(aux) else rownames(aux)
    cn = if(is.null(colnames(aux))) 1:ncol(aux) else colnames(aux)
    Dimnames=list(rn, cn, nm)
    if(!is.null(stackDimName)){
      if(is.character(stackDimName)){
        names(Dimnames) = c("loc","var", stackDimName)
        nmd = names(dimnames.data.frame(aux))
        tk = !is.null(nmd)
        if(any(tk))
          names(Dimnames)[tk] <- nmd[tk]
      }
    }
  } 
  uncolname = function(x){
    colnames(x) = NULL
    x
  }
  for(i in 2:length(x)){
    aux = cbind(aux,uncolname(x[[i]]))
  }
  aux = data.frame(aux)
  res = gsi.dataframestack(aux, stackDim=stackDimName, Dimnames=Dimnames,...)

  return(res)
}


#' @describeIn DataFrameStack.data.frame create a DataFrameStack from an array
#' @method as.DataFrameStack list
#' @export
as.DataFrameStack.list <- DataFrameStack.list


#' Convert a stacked data frame into an array
#' 
#' Recast a [DataFrameStack()] into a named array.
#'
#' @param x input [DataFrameStack()]
#' @param ... generic consistency
#'
#' @return the data recasted as an array with appropriate names.
#' @export
#'
#' @examples
#' ll = lapply(1:3, function(i) as.data.frame(matrix(1:10+(i-1)*10, ncol=2)))
#' dfs = DataFrameStack(ll, stackDimName="rep")
#' as.array(dfs)
#' 
as.array.DataFrameStack = function(x, ...){
  aux = x
  x = as.matrix(x)
  dim(x) = sapply(attr(aux,"Dimnames"), length)
  dimnames(x) = dimnames(aux)
  class(x) = "array"
  return(x)
}

#' Convert a stacked data frame into a list of data.frames
#' 
#' Recast a [DataFrameStack()] into a list of data.frames
#'
#' @param x input [DataFrameStack()]
#' @param ... generic consistency
#'
#' @return the data recasted as list of data.frames
#' @export
#'
#' @examples
#' ar = array(1:30, dim = c(5,2,3), dimnames=list(obs=1:5, vars=c("A","B"), rep=1:3))
#' dfs = DataFrameStack(ar, stackDim="rep")
#' as.list(dfs)
as.list.DataFrameStack = function(x, ...){
  res = lapply(dimnames(x)[attr(x, "stackDim")][[1]], function(i)
    getStackElement(x,i)
  )
  nm = dimnames(x)[stackDim(x)]
  if(is.null(nm)) return(res)
  names(res) = nm[[1]]
  return(res)
}

#' Return the dimnames of a DataFrameStack
#' 
#' Return the dimnames of a [DataFrameStack()], i.e. the three dimensions
#' 
#' @param x a [DataFrameStack()] object
#' 
#' @return a list (possibly named) with the element names of each of the three dimensions 
#' 
#' @export
dimnames.DataFrameStack = function(x) attr(x, "Dimnames")


#' @describeIn dimnames.DataFrameStack dimnames method for all Spatial*DataFrame objects of 
#' package `sp` which data slot contains a [DataFrameStack()]
#' @method dimnames Spatial
setMethod("dimnames", signature="Spatial", definition = function(x) dimnames(x@data) )




#' Get/set name/index of (non)stacking dimensions
#' 
#' Return (or set) the name or index of either the stacking dimension, or else of the 
#' non-stacking dimension (typically, the dimension runing through the variables)
#'
#' @param x a [DataFrameStack()] object, (only for `stackDim` it can also be a `Spatial` 
#' object which `data` slot is a `DataFrameStack`)
#' @param ... extra arguments for generic functionality
#'
#' @return the index or the name of the asked dimension.
#' @export
#'
#' @examples
#' ar = array(1:30, dim = c(5,2,3), dimnames=list(obs=1:5, vars=c("A","B"), rep=1:3))
#' dfs = DataFrameStack(ar, stackDim="rep")
#' dfs
#' stackDim(dfs)
#' noStackDim(dfs)
#' getStackElement(dfs, 1)
#' stackDim(dfs) <- "vars"
#' getStackElement(dfs, 1)
stackDim <- function(x,...) UseMethod("stackDim",x)


# @describeIn stackDim Get/set name/index of (non)stacking dimensions
#' @method stackDim default
#' @export
stackDim.default = function(x, ...) attr(x, "stackDim")


#' @describeIn stackDim Get/set name/index of (non)stacking dimensions
#' @method stackDim DataFrameStack
#' @export
stackDim.DataFrameStack = function(x, ...) attr(x, "stackDim")

setGeneric("stackDim", def=function(x, ...) attr(x, "stackDim") )

#' Get name/index of the stacking dimension of a Spatial object
#' 
#' Get name/index of the stacking dimension of a Spatial object which 
#' data slot is of class [DataFrameStack()]
#' 
#' @seealso [stackDim()]
#'  
#' @method stackDim Spatial
#' @param x a `Spatial` object which `data` slot is a `DataFrameStack`
#' @export
#' @return see [stackDim()] for details
#' @importFrom sp coordinates CRS bbox
setMethod("stackDim", signature="Spatial", definition = function(x) stackDim(x@data) )


#' @describeIn stackDim Get/set name/index of (non)stacking dimensions
#' @export
noStackDim <- function(x,...) UseMethod("noStackDim",x)

#' @describeIn stackDim Get/set name/index of (non)stacking dimensions
#' @method noStackDim default
#' @export
noStackDim.default = function(x, ...){
  if(stackDim(x)==2)
    return(3)
  if(stackDim(x)==names(dimnames(x))[2])
    return(names(dimnames(x))[3])
  if(stackDim(x)==3)
    return(2)
  if(stackDim(x)==names(dimnames(x))[3])
    return(names(dimnames(x))[2])
  warning("noStackDim: stackDim not found or not understood!")
  return(NA)
}


#' @describeIn stackDim Get/set name/index of (non)stacking dimensions
#' @param value the name or the index to be considered as stacking dimension
#' @export
`stackDim<-` <- function(x, value) UseMethod("stackDim<-", x)


#' @describeIn stackDim Get/set name/index of (non)stacking dimensions
#' @method stackDim<- default
#' @export
`stackDim<-.default` = function(x, value){
  attr(x, "stackDim") <- value
  return(x)
} 

#' Extract rows of a DataFrameStack
#' 
#' Extract rows (and columns if you know what you are doing) from a stacked data frame 
#'
#' @param x [DataFrameStack()] object
#' @param i row indices, names or logical vector of appropriate length (i.e. typically locations, observations, etc)
#' @param j column indices, names or logical vector of appropriate length. DO NOT USE if you are 
#' not sure what you are doing. The result will be a conventional data.frame, probably with the 
#' stacking structure destroyed.
#' @param ... generic parameters, ignored.
#' @param drop logical, if selection results in one single row or column selected, return a vector?
#'
#' @return the DataFrameStack or data.frame of the selected rows and columns.
#' @export
#'
#' @examples
#' ar = array(1:30, dim = c(5,2,3), dimnames=list(obs=1:5, vars=c("A","B"), rep=1:3))
#' dfs = DataFrameStack(ar, stackDim="rep")
#' dfs[1:2,]
#' stackDim(dfs[1:2,])
`[.DataFrameStack` <- function(x,i=NULL,j=NULL,..., drop=FALSE){
  erg = x
  aux = attr(x,"Dimnames")
  rnms = aux[[1]]
  class(erg) = "data.frame"
  if(!is.null(i)){
    erg = erg[i,, drop=drop]
    if(is.numeric(i) | is.logical(i)){
      aux[[1]] = rnms[i]
    }else if(is.character(i)){
      aux[[1]] = i
    }
    if(is.null(j)) return(gsi.dataframestack(erg, stackDim=stackDim(x), Dimnames=aux))
  }
  return(erg[,j, drop=drop])
}


#' Set or get the i-th data frame of a data.frame stack
#' 
#' Set or get one element of the [DataFrameStack()] 
#'
#' @param x container data, typically a [DataFrameStack()], but it can also be certain [sp::Spatial()]
#' object derivates of it
#' @param i index (or name) of the element of the stack to extract or replace 
#' @param ... extra arguments for generic functionality
#'
#' @return For the getters, the result is the data.frame of the stack asked for. For the setters
#' the result is the original DataFrameStack with the corresponding element replaced. Spatial methods
#' return the corresponding spatial object, ie. the spatial information of the stack is transferred to the
#' extracted element.
#' @export
#' @aliases setStackElement
#' @examples
#' ar = array(1:30, dim = c(5,2,3), dimnames=list(obs=1:5, vars=c("A","B"), rep=1:3))
#' dfs = DataFrameStack(ar, stackDim="rep")
#' dfs
#' stackDim(dfs)
#' getStackElement(dfs, 1)
getStackElement <- function(x,i,...) UseMethod("getStackElement",x)


#' @describeIn getStackElement Set or get one element of the [DataFrameStack()] 
#' @method getStackElement default
#' @export
getStackElement.default = function(x, i, ...){
  if(is(x, "SpatialPixelsDataFrame"))return(getStackElement_SpatialPixelsDataFrame(x=x, i=i, ...))
  if(is(x, "SpatialGridDataFrame"))return(getStackElement_SpatialGridDataFrame(x=x, i=i, ...))
  if(is(x, "SpatialPointsDataFrame"))return(getStackElement_SpatialPointsDataFrame(x=x, i=i, ...))
  x[i, drop=F]
} 


#' @describeIn getStackElement Set or get one element of the [DataFrameStack()] 
#' @method getStackElement list
#' @export
getStackElement.list = function(x, i, ...) x[[i]]



#' @describeIn getStackElement Set or get one element of the [DataFrameStack()] 
#' @method getStackElement DataFrameStack
#' @param MARGIN which dimension is the stacking dimension? you seldom want to touch this!!
#' @export
getStackElement.DataFrameStack = function(x, i, MARGIN=stackDim(x), ...){
  dm = sapply(attr(x,"Dimnames")[-1], length)
  l = prod(dm)
  idx = matrix(1:l, ncol=dm[1], nrow=dm[2], byrow=TRUE)
  colnames(idx) = attr(x,"Dimnames")[[2]]
  rownames(idx) = attr(x,"Dimnames")[[3]]
  if((stackDim(x)==2) | (stackDim(x)==names(dimnames(x))[2])){
    erg = data.frame(x[,idx[,i], drop=F])
    colnames(erg) = attr(x,"Dimnames")[[3]]
  }else{
    erg = data.frame(x[,idx[i,], drop=F])
    colnames(erg) = attr(x,"Dimnames")[[2]]
  } 
  rownames(erg) = attr(x,"Dimnames")[[1]]
  
  class(erg) = setdiff(class(x), "DataFrameStack")
  return(erg)
}


## getStackElement S4 methods:
# @rdname stackElementsS4 get/setStackElement S4 methods
# 
# These functions provide S4 methods for the [getStackElement()]/`SetStackElement()` 
# generic functionality
# @export
# setGeneric("getStackElement", def=function(x, ...) standardGeneric("getStackElement"))

getStackElement_SpatialGridDataFrame = function(x,i,MARGIN=stackDim(x), ...){
  sp::SpatialGridDataFrame(grid=sp::getGridTopology(x), data = getStackElement(x@data, i),
                       proj4string = sp::CRS(sp::proj4string(x)))
}
# @rdname stackElementsS4 Get method for SpatialGridDataFrame object
# @method getStackElement SpatialGridDataFrame
# @export
#setMethod("getStackElement", signature="SpatialGridDataFrame", def=getStackElement_SpatialGridDataFrame)


getStackElement_SpatialPointsDataFrame = function(x,i,MARGIN=stackDim(x), ...){
  sp::SpatialPointsDataFrame(coords=sp::SpatialPoints(sp::coordinates(x)), 
                             data = getStackElement(x@data, i),
                       proj4string = sp::CRS(sp::proj4string(x)), 
                       bbox = sp::bbox(x))
}
# @rdname stackElementsS4 Get method for SpatialPointsDataFrame object
# @method getStackElement SpatialPointsDataFrame
# @export
#setMethod("getStackElement", signature="SpatialPointsDataFrame", def=getStackElement_SpatialPointsDataFrame)


getStackElement_SpatialPixelsDataFrame = function(x,i,MARGIN=stackDim(x), ...){
  sp::SpatialPixelsDataFrame(points=sp::SpatialPoints(sp::coordinates(x)), 
                             data = getStackElement(x@data, i),
             proj4string = sp::CRS(sp::proj4string(x)), 
             grid = sp::getGridTopology(x))
}
# @rdname stackElementsS4 Get method for SpatialPixelsDataFrame object
# @method getStackElement SpatialPixelsDataFrame
# @export
# setMethod("getStackElement", signature="SpatialPixelsDataFrame", def=getStackElement_SpatialPixelsDataFrame)



# @describeIn getStackElement Set one element of the [DataFrameStack()] 
#' @export
setStackElement <- function(x,i,...) UseMethod("setStackElement",x)


#' @describeIn getStackElement Set or get one element of the [DataFrameStack()] 
#' @method setStackElement default
#' @param value for the setting operation, the new data.frame to replace the selected one; note
#' that the compatibility of the dimensions of `value` is only checked for setStackElement.DataFrameStack
#' and its Spatial derivates; for other methods `setStackElement` can break the consistency of the  
#' stack!
#' @export
setStackElement.default = function(x, i, value, ...){
  if(is(x, "SpatialPixelsDataFrame"))return(setStackElement_SpatialPixelsDataFrame(x=x, i=i, value=value, ...))
  if(is(x, "SpatialGridDataFrame"))return(setStackElement_SpatialGridDataFrame(x=x, i=i, value=value, ...))
  if(is(x, "SpatialPointsDataFrame"))return(setStackElement_SpatialPointsDataFrame(x=x, i=i, value=value, ...))
  x[i] <- value
  return(x)
} 

#' @describeIn getStackElement Set  one element of the [DataFrameStack()] in data.frame form 
#' @method setStackElement data.frame
#' @export
setStackElement.data.frame = function(x, i, value, ...){
  if(is(x, "SpatialPixelsDataFrame"))return(setStackElement_SpatialPixelsDataFrame(x=x, i=i, value=value,...))
  if(is(x, "SpatialGridDataFrame"))return(setStackElement_SpatialGridDataFrame(x=x, i=i, value=value,...))
  if(is(x, "SpatialPointsDataFrame"))return(setStackElement_SpatialPointsDataFrame(x=x, i=i, value=value,...))
  
  x[,i] <- value
  return(x)
} 

#' @describeIn getStackElement Set get one element of a [DataFrameStack()] in list form
#' @method setStackElement list
#' @export
setStackElement.list = function(x, i, value, ...){
  x[[i]] <- value
  return(x)
}

#' @describeIn getStackElement Set one element of the [DataFrameStack()] 
#' @method setStackElement DataFrameStack
#' @export
setStackElement.DataFrameStack = function(x,i, value, MARGIN=stackDim(x), ...){
  dm = sapply(attr(x,"Dimnames")[-1], length)
  l = prod(dm)
  idx = matrix(1:l, ncol=dm[1], nrow=dm[2], byrow=TRUE)
  aux = getStackElement(x,i, MARGIN=MARGIN)
  colnames(idx) = attr(x,"Dimnames")[[2]]
  rownames(idx) = attr(x,"Dimnames")[[3]]
  # do some checks
  if(length(dim(value))==0)  dim(value) = c(length(value),1)
  if(any(dim(value)!=dim(aux))) stop("setStackElement.DataFrameStack: provided value does not have the appropriate dimensions")
  if(any(colnames(value)!=colnames(aux))) warnings("setStackElement.DataFrameStack: colnames of provided value do not fully correspond with stack varnames")
  # replace and return
  if((MARGIN==2) | (MARGIN==names(dimnames(x))[2])){
    x[,idx[,i]] <- value
  }else{
    x[,idx[i,]] <- value
  } 
  return(x)
}

## getStackElement S4 methods:
# @rdname stackElementsS4 Set methods for S4 classes
# @export
#setGeneric("setStackElement", def=function(x, ...) standardGeneric("setStackElement"))


setStackElement_SpatialGridDataFrame = function(x,i, value, MARGIN=stackDim(x), ...){
  sp::SpatialGridDataFrame(grid=sp::getGridTopology(x), data = setStackElement(x@data, i, value),
                       proj4string = sp::CRS(sp::proj4string(x)))
}
# @rdname stackElementsS4 Set method for SpatialGridDataFrame object
# @method setStackElement SpatialGridDataFrame
# @export
#setMethod("setStackElement", signature="SpatialGridDataFrame", def=setStackElement_SpatialGridDataFrame)


setStackElement_SpatialPointsDataFrame = function(x,i, value, MARGIN=stackDim(x), ...){
  sp::SpatialPointsDataFrame(coords=sp::coordinates(x), data = setStackElement(x@data, i, value),
                         proj4string = sp::CRS(sp::proj4string(x)), bbox = bbox(x))
}
# @rdname stackElementsS4 Set method for SpatialPointsDataFrame object
# @method setStackElement SpatialPointsDataFrame
# @export
#setMethod("setStackElement", signature ="SpatialPointsDataFrame", def=setStackElement_SpatialPointsDataFrame)


setStackElement_SpatialPixelsDataFrame = function(x,i, value, MARGIN=stackDim(x), ...){
  sp::SpatialPixelsDataFrame(points=sp::coordinates(x), data = setStackElement(x@data, i, value),
                       proj4string = sp::CRS(sp::proj4string(x)), grid = sp::getGridTopology(x))
}
# @rdname stackElementsS4 Set method for SpatialPixelsDataFrame object
# @method setStackElement SpatialPixelsDataFrame
# @export
#setMethod("setStackElement", signature="SpatialPixelsDataFrame", def=setStackElement_SpatialPixelsDataFrame)



### tests -----------
# if(!exists("do.test")) do.test=FALSE
# if(do.test){
#   
#   mm = array(rnorm(60), dim=c(3,5,4), dimnames=list(loc=letters[1:3], sim=1:5, var=LETTERS[1:4]))
#   mdt = data.frame(matrix(mm, nrow=3, ncol=20))
#   dfs = DataFrameStack(mdt, dimnames=list(loc=letters[1:3], sim=1:5, var=LETTERS[1:4])) 
#   dimnames(dfs)
#   as.array(dfs)
#   stackDim(dfs)
#   a = SpatialPointsDataFrame(coords=expand.grid(x=-1:1, y=0), data = dfs)
#   stackDim(a)
#   dimnames(a)
#   as.list(dfs)
#   plot(a)
#   spplot(a)
#   a@data
#   class(a@data)
# }
  



  
#' Plot a swarm of calculated output through a DataFrameStack
#' 
#' Take a [DataFrameStack()] and apply a certain plotting function to each elements of the stack.
#' The result (typically a curve per each stack element), may then be plotted all together
#'
#' @param X a [DataFrameStack()] or anaologous object
#' @param MARGIN which dimension defines the stack? it has a good default! change only if you
#' know what you do 
#' @param PLOTFUN the elemental calculating function; this must take as input one element
#' of the stack and return as output the (x,y)-coordinates of the calculated curve for that
#' element, in a list of two elements
#' @param ... further parameters to `PLOTFUN`
#' @param .plotargs either a logical, or else a list of graphical arguments to pass
#' to [plot.swarmPlot()]; if `.plotargs=FALSE` no plot is produced; if `.plotargs=TRUE` 
#' a plot with default values is produced; 
#' @param .parallelBackend NA or a parallelization strategy; currently unstable for certain 
#' operations and platforms. 
#'
#' @return Invisibly, this function returns a list of the evaluation of `PLOTFUN` on 
#' each element of the stack. If `.plotargs` other than `FALSE`, then the function calls 
#' [plot.swarmPlot()] to produce a plot.
#' @export
#' @import foreach
#' @importFrom graphics lines.default points.default
#'
#' @examples
#' dm = list(point=1:100, var=LETTERS[1:2], rep=paste("r",1:5, sep=""))
#' ar = array(rnorm(1000), dim=c(100,2,5), dimnames = dm)
#' dfs = DataFrameStack(ar, stackDim="rep")
#' swarmPlot(dfs, PLOTFUN=function(x) density(x[,1]))
swarmPlot = function(X, MARGIN=stackDim(X), PLOTFUN,..., 
                     .plotargs=list(type="l"),
                     .parallelBackend=NA){
  PLTFN = function(i,...) PLOTFUN(getStackElement(X,i=i),...)
  if(is.na(.parallelBackend)){
    idx = 1:length(dimnames(X)[[MARGIN]])
    resL = lapply(idx, FUN=PLTFN, ...)
  }else{
    requireNamespace("foreach", quietly = TRUE)
    message("swarmPlot: it is your responsibility to open, provide and close the right .parallelBackend!")
    idx = 1:length(dimnames(X)[[MARGIN]])
    resL = foreach(idx=idx) %dopar% PLTFN(idx)
  }  
  if(is.logical(.plotargs))
    if(!.plotargs){
      class(resL) = "swarmPlot"
      return(invisible(resL))
    }else{
      .plotargs = list()
    }
  class(resL) = "swarmPlot"
  .plotargs$x = resL
  do.call(plot.swarmPlot, .plotargs)
  invisible(resL)
}



#' Plotting method for swarmPlot objects
#'
#' Produces a plot for objects obtained by means of [swarmPlot()]. These are lists of two-element
#' list, respectively giving the (x,y) coordinates of a series of points defining a curve. They
#' are typically obtained out of a [DataFrameStack()], and each curve represents a a certain relevant
#' summary of one of the elements of the stack. All parameters of [graphics::plot.default()] 
#' are available, plus the following
#'
#' @param x swarmPlot object
#' @param type see [graphics::plot.default()]
#' @param col color to be used for all curves and points
#' @param xlim limits of the X axis, see [graphics::plot.default()]
#' @param ylim limits of the Y axis, see [graphics::plot.default()]
#' @param xlab label for the X axis, see [graphics::plot.default()]
#' @param ylab label for the Y axis, see [graphics::plot.default()]
#' @param main main plot title, see [graphics::plot.default()]
#' @param asp plot aspect ratio, see [graphics::plot.default()]
#' @param sub plot subtitle, see  [graphics::plot.default()]
#' @param log log scale for any axis? see [graphics::plot.default()]
#' @param ann plot annotation indications, see [graphics::plot.default()]
#' @param axes should axes be set? see [graphics::plot.default()]
#' @param frame.plot should the plot be framed? see [graphics::plot.default()]
#' @param bg fill color for the data points if certain `pch` symbols are used
#' @param pch symbol for the data points, see [graphics::points()]
#' @param cex symbol size for the data points
#' @param lty line style for the curves
#' @param lwd line thickness for the curves
#' @param add logical, add results to an existing plot?
#' @param ... further parameters to the plotting function
#'
#' @return Nothing. The function is called to make a plot.
#' @export
#' @importFrom graphics plot lines points
#'
#' @examples
#' dm = list(point=1:100, var=LETTERS[1:2], rep=paste("r",1:5, sep=""))
#' ar = array(rnorm(1000), dim=c(100,2,5), dimnames = dm)
#' dfs = DataFrameStack(ar, stackDim="rep")
#' rs = swarmPlot(dfs, PLOTFUN=function(x) density(x[,1]), .plotargs=FALSE)
#' plot(rs, col="yellow", lwd=3)
plot.swarmPlot = function(x, type="l", col="grey", xlim=NULL, ylim=NULL, xlab="x", ylab="y", 
                          main=NULL, asp=NA, sub=NULL, log="", ann = par("ann"), axes = TRUE, 
                          frame.plot = axes,
                          bg=NA, pch=1, cex=1, lty=1, lwd=1, add=FALSE, ...){
  resL = x
  idx = 1:length(x)
  if(is.null(xlim))
    xlim = range(sapply(resL, function(x) x[[1]]))
  if(is.null(ylim))
    ylim = range(sapply(resL, function(x) x[[2]]))
  if(!add){
    x = resL[[1]][[1]]
    y = resL[[1]][[2]]
    plot(x=x, y=y, type=type, col=col, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, 
         asp=asp, sub=sub, log=log, axes=axes, frame.plot=frame.plot, bg=bg, pch=pch, cex=cex,
         lty=lty, lwd=lwd, ...)
    idx = idx[-1]
  }
  for(i in idx){
    x = resL[[i]][[1]]
    y = resL[[i]][[2]]
    if(type %in% c("l", "b", "o"))
      lines.default(x=x, y=y,col=col, lwd=lwd, lty=lty)
    if(type %in% c("p", "b", "o"))
      points.default(x=x, y=y,col=col, pch=pch, cex=cex, bg=bg)
  }
}

#' Apply Functions Over Array or DataFrameStack Margins
#' 
#' Returns a vector or array or list of values obtained by 
#' applying a function to the margins of an array or matrix.
#' Method `gmApply.default()` is just a wrapper on [base::apply()].
#' Method `gmApply()` reimplements the functionality
#' with future access to parallel computing and appropriate default 
#' values for the MARGIN. ALWAYS use named arguments here!
#' 
#' @param X a [DataFrameStack()] object (see [base::apply()] for other options)
#' @param ... further arguments to `FUN`
#' 
#' @return In principle, if `MARGIN==stackDim(X)` (the default), the oputput is a list
#' with the result of using `FUN` on each element of the stack. If `FUN` returns a 
#' matrix or a data.frame assimilable to one element of the stack, a transformation of
#' this output to a DataFrameStack is attempted. 
#' 
#' For `X` non-DataFrameStack or `MARGIN!=stackDim(X)` see [base::apply()].
#' @export
gmApply <- function(X, ...) UseMethod("gmApply", X)


#' @describeIn gmApply wrapper around [base::apply()]
#' @param MARGIN a name or an index of the dimension along which should 
#' the calculations be done; defaults to the stacking dimension of the 
#' [DataFrameStack()], i.e. to the output of `stackDim(X)`
#' @param FUN function to apply; the default behaviour being that this function
#' is applied to each element of the stack `X`
#' @param ... further arguments to `FUN`
#'
#' @export
#' @method gmApply default
gmApply.default <- function(X, MARGIN, FUN, ...) base::apply(X, MARGIN, FUN, ...)


#' @describeIn gmApply Apply Functions Over DataFrameStack Margins
#' @param .parallel currently ignored
#' @method gmApply DataFrameStack
#' @export
#' @examples
#' dm = list(point=1:100, var=LETTERS[1:2], rep=paste("r",1:5, sep=""))
#' ar = array(rnorm(1000), dim=c(100,2,5), dimnames = dm)
#' dfs = DataFrameStack(ar, stackDim="rep")
#' gmApply(dfs, FUN=colMeans)
#' rs = gmApply(dfs, FUN=function(x) x+1)
#' class(rs)
#' getStackElement(rs,1)
#' getStackElement(dfs,1)
gmApply.DataFrameStack = function(X, MARGIN=stackDim(X), FUN, ..., .parallel=FALSE){
  # manage the general case of apply
  if(any(MARGIN!=stackDim(X))){
    return(apply(as.array(X), MARGIN=MARGIN, FUN=FUN, ...)) 
  }
  # manage the special case of fun across the stacking dimension, much cheaper:
  FN = function(i) FUN(getStackElement(X,i=i, MARGIN=MARGIN),...)
  # do the calculations
  if(.parallel){
    stop(".parallel not yet implemented")
  }else{
    idx = 1:length(dimnames(X)[[MARGIN]])
    res = lapply(idx, FUN=FN)
  }  
  if(is.numeric(MARGIN))
    MARGIN = names(attr(X, "Dimnames"))[MARGIN]
  # redimensionalise and cast to DataFrameStack
  if(is.character(MARGIN)){
    out = DataFrameStack(res, stackDimName = MARGIN, ...)
  }else{
    out = DataFrameStack(res, ...)
  }
  # pass all attributes:
  #attr(out,"Dimnames") = list(dimnames(res[[1]], )
  out = gsi.matchDimnames(X,out)
  attrb = names(attributes(X))
  attrb = setdiff(attrb, c("names", "col.names", "row.names", "dimnames", "Dimnames", "dim", "class"))
  if(length(attrb)>0)
    for(atr in attrb) attr(out,atr) = attr(X, atr)
  return(out)
}

gsi.matchDimnames = function(X, out){
  dimnO = dimnames(out)
  dmX = sapply(dimnames(X), length)
  dmO = sapply(dimnO, length)
  oo = outer(dmX, dmO, "==")
  if(any(oo)){
    ii = which(colSums(oo)==1)
    jj = unlist(apply(oo,2,which))
    for(i in ii){
      j = which(oo[,i])
      dimnO[[i]]=dimnames(X)[[j]]
    }
    names(dimnO)[ii] <- names(dimnames(X))[jj]
    if(any(is.na(names(dimnO))))
      names(dimnO)[is.na(names(dimnO))] = paste("d",which(is.na(names(dimnO))), sep="")
    attr(out, "Dimnames") = dimnO
  }
  return(out)
}


### tests -----------
# if(!exists("do.test")) do.test=FALSE
# if(do.test){
#   
#   exampleFun = function(x){
#     print(x[,1])
#     density(x[,1], n=256)
#   } 
#   swarmPlot(dfs,PLOTFUN=exampleFun, .parallel = F)
#   
#   exampleFun = function(x){
#     erg = as.matrix(x)+10*rep(1:4, each=3)
#     dimnames(erg) = dimnames(x)
#     return(erg)
#   } 
#   aux = gmApply(dfs, FUN=exampleFun)
#   
#   library(gstat)
#   library(sp)
#   data("jura", package = "gstat")
#   coords = SpatialPoints(jura.pred[,1:2])
#   plot(coords)
#   
# }
