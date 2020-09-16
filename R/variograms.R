#### theoretical variogram ---------------------

#' gmGeostats Variogram models
#' set up a D-variate variogram models 
#'
#' @param type constant indicating the model of correlation function
#' @param nugget (DxD)-matrix for the nugget effect
#' @param sill (DxD)-matrix for the partial sills of the correlation function
#' @param anisRanges 2x2 or 3x3 matrix of ranges (see details)
#' @param extraPar for certain correlation functions, extra parameters (smoothness, period, etc)
#'
#' @return an object of class "gmCgram" containing the linear model of corregionalization
#' of the nugget and the structure given. 
#' @details The argument `type` must be an integer indicating the model to be used. 
#' Some constants are available to make reading code more understandable. That is, you can 
#' either write `1`, `vg.Sph` or `vg.Spherical`, they will all work and produce
#' a spherical model. The same applies for the following models: `vg.Gauss=0`; 
#' `vg.Exp=vg.Exponential=2`. These constants are available after calling
#' `data("variogramModels")`. 
#' No other model is currently available, but this data object will be 
#' regularly updated.
#' The constant vector `gsi.validModels` contains all currently valid models.
#' 
#' Argument `anisRange` expects a matrix $M$ such that 
#' \deqn{
#' h^2 = (\mathbf{x}_i-\mathbf{x}_j)\cdot M^{-1}\cdot (\mathbf{x}_i-\mathbf{x}_j)^t
#' }
#' is the (square of) the lag distance to be fed into the correlation function.
#' @family gmCgram 
#' @export
#' @aliases vg.Exp vg.exp vg.Exponential vg.Gau vg.gauss 
#' vg.Gauss vg.Sph vg.sph vg.Spherical gsi.validModels
#'
#' @examples
#' utils::data("variogramModels") # shortcut for all model constants
#' v1 = setCgram(type=vg.Gau, sill=diag(2), anisRanges = 3*diag(c(3,1)))
#' v2 = setCgram(type=vg.Exp, sill=0.3*diag(2), anisRanges = 0.5*diag(2))
#' vm = v1+v2
#' plot(vm)
setCgram = function(type, nugget=sill*0, sill, anisRanges, extraPar=0){
  utils::data("variogramModels")
  stopifnot(all(dim(sill)==dim(nugget)), 
            ncol(sill)==nrow(sill),
            type %in% gsi.validModels)
  dim(sill) = c(1,dim(sill))
  dim(anisRanges) = c(1,dim(anisRanges))
  vgout <- list(type=type,
                data=extraPar,
                nugget=nugget,
                sill=sill, #(nstru, nvar, nvar)
                M=anisRanges     # these are "classical" ranges (i.e. distances)
  )
  class(vgout) = "gmCgram"
  return(vgout)
}



#' Subsetting of gmCgram variogram structures
#' 
#' Extraction or combination of nested structures of a gmCgram object
#'
#' @param x `gmCgram` variogram object
#' @param i indices of the structures that are desired to be kept (0=nugget) or removed (see details)
#' @param ... extra arguments for generic functionality
#'
#' @return a `gmCgram` variogram object with the desired structures only. 
#' @details This function can be used to: extract the nugget (i=0), extract some 
#' structures (i=indices of the structures, possibly including 0 for the nugget),
#' or filter some structures out (i=negative indices of the structures to remove; 
#' nugget will always removed in this case). If you want to extract "slots" or 
#' "elements" of the variogram, use the $-notation. If you want to extract variables of the
#' variogram matrices, use the `[`-notation. The contrary operation (adding structures together)
#' is obtained by summing (+) two `gmCgram` objects.
#' @export
#' @method [[ gmCgram
#' @family gmCgram functions
#' @examples
#' utils::data("variogramModels")
#' v1 = setCgram(type=vg.Gau, sill=diag(2), anisRanges = 3*diag(c(3,1)))
#' v2 = setCgram(type=vg.Exp, sill=0.3*diag(2), anisRanges = 0.5*diag(2))
#' vm = v1+v2
#' vm[[1]]
"[[.gmCgram"<- function(x,i,...){
  nG = dim(x$M)[3]
  nD = dim(x$nugget)[1]
  nSo = length(i)
  nSi = dim(x$M)[1]
  if(i==0){# case only nugget wanted
    out = with(x, list(type=type[i], 
                       data=data[i], 
                       nugget=nugget, 
                       sill=structure(0*sill[1,,], dim=c(nSo, nD, nD)), 
                       M=structure(M[1,,], dim=c(nSo, nG, nG))))
  }else if( (0 %in% i) & !any(i<0)){
    j = i[i>0] # case some structures wanted, including nugget
    out = with(x, list(type=type[j], 
                       data=data[j], 
                       nugget=nugget, 
                       sill=structure(sill[j,,], dim=c(nSo-1, nD, nD)), 
                       M=structure(M[j,,], dim=c(nSo-1, nG, nG))))
  }else if(all(i>0) | all(i<0)){# case some structured (un)wanted, nugget surely not wanted
    if(all(i<0)) nSo = nSi-nSo
    out = with(x, list(type=type[i], 
                       data=data[i], 
                       nugget=nugget*0, 
                       sill=structure(sill[i,,], dim=c(nSo, nD, nD)), 
                       M=structure(M[i,,], dim=c(nSo, nG, nG))))
  }else{# unsolvable case
    stop("index set i cannot merge 0 and negative numbers")
  }
  class(out) = class(x)
  return(out)
}


#' Combination of gmCgram variogram structures
#' 
#' combination of nested structures of a gmCgram object
#'
#' @param x `gmCgram` variogram object
#' @param y `gmCgram` variogram object
#' @export
#' @return The combined nested structures
#' @method + gmCgram
#' @examples
#' utils::data("variogramModels")
#' v1 = setCgram(type=vg.Gau, sill=diag(2), anisRanges = 3*diag(c(3,1)))
#' v2 = setCgram(type=vg.Exp, sill=0.3*diag(2), anisRanges = 0.5*diag(2))
#' vm = v1+v2
"+.gmCgram" <- function(x,y) {
  y = as.gmCgram(y)
  stopifnot(class(y)=="gmCgram", 
            dim(x$sill)[-1]==dim(y$sill)[-1],
            dim(x$M)[-1]==dim(y$M)[-1])
  myfun = function(A,B){
    D = dim(A)[2]
    nA = dim(A)[1]
    nB = dim(B)[1]
    dim(A) = c(nA, D^2)
    dim(B) = c(nB, D^2)
    out = rbind(A, B)
    dim(out) = c(nA+nB, D, D)
    return(out)
  }
  x$type = c(x$type, y$type)
  x$data= c(x$data, y$data)
  x$nugget = x$nugget + y$nugget
  x$sill = myfun(x$sill, y$sill)
  x$M = myfun(x$M, y$M)
  return(x)
}




#' Subsetting of gmCgram variogram structures
#' 
#' Extraction of some variables of a gmCgram object
#'
#' @param x \code{gmCgram} variogram object
#' @param i row-indices of the variables to be kept/removed
#' @param j column-indices of the variables to be kept/removed (if only \code{i}
#' is specified, \code{j} will be taken as equal to \code{i}!)
#' @param ...  extra arguments for generic functionality
#'
#' @return a \code{gmCgram} variogram object with the desired variables only. 
#' @details This function can be used to extract the model for a a subset of variables. 
#' If only \code{i} is specified, \code{j} will be taken as equal to \code{i}. 
#' If you want to select all \code{i}'s for certain \code{j}'s or vice versa, give
#' \code{i=1:dim(x$nugget)[1]} and \code{j=} your desired indices, respectively  
#' \code{j=1:dim(x$nugget)[2]} and \code{i=} your desired indices; replace \code{x} by the
#' object you are giving. If \code{i!=j}, the output will be a \code{c("gmXCgram","gmCgram")} 
#' object, otherwise it will be a regular class \code{"gmCgram"} object.
#' If you want to extract "slots" or 
#' "elements" of the variogram, use the $-notation. If you want to extract variables of the
#' variogram matrices, use the `[`-notation.
#' @export
#' @method [ gmCgram
#' @family gmCgram functions
#' @examples
#' utils::data("variogramModels")
#' v1 = setCgram(type=vg.Gau, sill=diag(2), anisRanges = 3*diag(c(3,1)))
#' v2 = setCgram(type=vg.Exp, sill=0.3*diag(2), anisRanges = 0.5*diag(2))
#' vm = v1+v2
#' vm[1,1]
"[.gmCgram"<- function(x,i,j=i,...){
  nDi = length(i)
  nDj = length(j)
  nS = dim(x$M)[1]
  out = with(x, list(type=type, 
                     data=data, 
                     nugget=structure(nugget[i,j, drop=F], dim=c(nDi, nDj)), 
                     sill=structure(sill[,i,j, drop=F], dim=c(nS, nDi, nDj)), 
                     M=M))
  class(out) = class(x)
  if(!all(i==j)) class(out) = unique(c("gmXCgram", class(out)))
  return(out)
}

#' Length, and number of columns or rows
#' 
#' Provide number of structures, and nr of variables of an LMC of class gmCgram 
#'
#' @param x gmCgram object
#'
#' @return \code{length} returns the number of structures (nugget not counted), while
#' \code{ncol} and \code{nrow} return these values for the nugget (assuming that they will
#' be also valid for the sill).
#' @export
#' @aliases ncol.gmCgram nrow.gmCgram
#' @family gmCgram functions
#'
#' @method length gmCgram
#' @examples
#' utils::data("variogramModels")
#' v1 = setCgram(type=vg.Gau, sill=diag(3)+0.5, anisRanges = 2*diag(c(3,0.5)))
#' v2 = setCgram(type=vg.Exp, sill=0.3*diag(3), anisRanges = 0.5*diag(2))
#' vm = v1+v2
#' length(vm)
#' ncol(vm)
#' nrow(vm)
length.gmCgram = function(x)  length(x$type)
if(!isGeneric("ncol")){
  ncol <- function(x) UseMethod("ncol",x)
  ncol.default <- base::ncol
}
if(!isGeneric("nrow")){
  nrow <- function(x) UseMethod("nrow",x)
  nrow.default <- base::nrow
}
ncol.gmCgram = function(x) ncol(x$nugget)
nrow.gmCgram = function(x) nrow(x$nugget)


#' Convert a gmCgram object to an (evaluable) function
#' 
#' Evaluate a gmCgram on some h values, or convert the gmCgram object into an evaluable function
#'
#' @param x a gmCgram object
#' @param ... extra arguments for generic functionality
#'
#' @return a \code{function} that can be evaluated normally, with an argument \code{X}
#' and possibly another argument \code{Y}; both must have the same number of columns 
#' than the geographic dimension of the variogram (i.e. \code{dim(x$M)[3]}).
#' @export
#' @method as.function gmCgram
#' @family gmCgram functions
#' @examples
#' utils::data("variogramModels")
#' v1 = setCgram(type=vg.Gau, sill=diag(2)+0.5, anisRanges = 2*diag(c(3,0.5)))
#' v2 = setCgram(type=vg.Exp, sill=0.3*diag(2), anisRanges = 0.5*diag(2))
#' vm = v1+v2
#' vgf = as.function(vm)
#' (h = rbind(c(0,1), c(0,0), c(1,1)))
#' vgf(h)
#' predict(vm, h)
as.function.gmCgram = function(x,...){
    f <- function(X,Y=X){
      if(is(X,"Spatial")) X = sp::coordinates(X)
      X = as.matrix(X)
      if(is(Y,"Spatial")) Y = sp::coordinates(Y)
      Y = as.matrix(Y)
      stopifnot(ncol(X)==ncol(Y), ncol(X)==dim(x$M)[3])
      ijEqual = ifelse(nrow(X)==nrow(Y), all(X==Y), FALSE)
      o = gsi.calcCgram(X,Y,x,ijEqual)
      return(o)
    }
  return(f)
}

#' @describeIn as.function.gmCgram predict a gmCgram object on some lag vector coordinates
#' @param object gmCgram object
#' @param newdata matrix, data.frame or Spatial object containing coordinates
#' @method predict gmCgram
#' @export
predict.gmCgram = function(object, newdata, ...){
  as.function(object)(X=newdata)
}


#' Convert theoretical structural functions to gmCgram format 
#' 
#' Convert covariance function or variogram models to the format gmCgram 
#' of package gmGeostats
#' 
#' @param m model to be converted
#' @param ... further parameters
#'
#' @return the covariance/variogram model, recasted to class \code{gmCgram}. 
#' This is a generic function. Methods exist for objects of class
#' \code{LMCAnisCompo} (for compositional data) and \code{variogramModelList}
#' (as provided by package \code{gstat}).  
#' @export
#' @family gmCgram functions
as.gmCgram <- function(m, ...)  UseMethod("as.gmCgram",m) 

#' @describeIn as.gmCgram Convert theoretical structural functions to gmCgram format 
#' @method as.gmCgram default
#' @export
as.gmCgram.default <- function(m,...) m 



#' Draw cuves for covariance/variogram models
#' 
#' Represent a gmCgram object as a matrix of lines in several plots
#'
#' @param x object to draw, of class gmCgram // curently only valid for symmetric functions
#' @param xlim.up range of lag values to use in plots of the upper triangle
#' @param xlim.lo range of lag values to use in plots of the lower triangle
#' @param vdir.up geograohic directions to represent in the upper triangle 
#' @param vdir.lo geograohic directions to represent in the lower triangle
#' @param xlength number of discretization points to use for the curves (defaults to 200)
#' @param varnames string vector, variable names to use in the labelling of axes
#' @param add logical, should a new plot be created or stuff be added to an existing one?
#' @param commonAxis logical, is a common Y axis for all plots in a row desired?
#' @param cov logical, should the covariance function (=TRUE) or the variogram (=FALSE) be plotted?   
#' @param closeplot logical, should the plot be left open (FALSE) for further changes, or be frozen (TRUE)? 
#' defaults to TRUE
#' @param ... further graphical parameters for the plotting function
#'
#' @return This function is called for its side effect of producing a plot: the plot will be
#' open to further changes if you provide `closeplot=FALSE`. Additionally, the function
#' invisibly returns the graphical parameters that were active before starting the plot. Hence, 
#' if you want to freeze a plot and not add anymore to it, you can do \code{par(plot(x, closeplot=FALSE, ...))},
#' or \code{plot(x, closeplot=TRUE, ...)}.
#' If you want to further add stuff to it, better just call \code{plot(x, closeplot=FALSE,...)}. The difference
#' is only relevant when working with the screen graphical device. 
#' @export 
#' @method plot gmCgram
#' @family gmCgram functions
#' @examples
#' utils::data("variogramModels")
#' v1 = setCgram(type=vg.Gau, sill=diag(3)-0.5, anisRanges = 2*diag(c(3,0.5)))
#' v2 = setCgram(type=vg.Exp, sill=0.3*diag(3), anisRanges = 0.5*diag(2))
#' vm = v1+v2
#' plot(vm)
#' plot(vm, cov=FALSE)
plot.gmCgram = function(x, xlim.up=NULL, xlim.lo=NULL, vdir.up= NULL, vdir.lo= NULL, xlength=200, varnames = colnames(x$nugget), 
                        add=FALSE, commonAxis=TRUE, cov =TRUE, closeplot=TRUE, ...){
  Dg = dim(x$M)[2]
  Ns = dim(x$M)[1]
  Dv = dim(x$nugget)[1]
  if(is.null(varnames)) varnames = paste("v", 1:Dv, sep="")
  if(is.null(vdir.up) & is.null(vdir.lo)){
    vdir.lo = rep(0, Dg)
    vdir.lo[1] = 1
    dim(vdir.lo) = c(1,Dg)
    if(!is.isotropic(x)){
      aux = diag(Dg)
      if(Dg==3){
        vdir.lo = aux[3,]
        vdir.up = aux[-3,]
      }else
        vdir.lo = aux
    }
  }
  if(!is.null(vdir.up)) vdir.up = compositions::oneOrDataset(compositions::normalize(vdir.up))
  if(!is.null(vdir.lo)) vdir.lo = compositions::oneOrDataset(compositions::normalize(vdir.lo))
  fk = c(sqrt(3),1,3)[x$type+1]*1.25  # range to effective range factor expansion
  if(is.null(xlim.up)){
    if(add){
      par(mfg=c(1,2))
      xlim.up = par()$usr[1:2]
    }else if(!is.null(vdir.up)){
      maxdist = sapply(1:Ns, function(i) max(sapply(1:nrow(vdir.up), function(j) vdir.up[j,]%*%x$M[i,,]%*%vdir.up[j,])))
      # here we must compute the max M projected on vd.up
      xlim.up = c(0, max(fk*maxdist))
    }
  }
  if(is.null(xlim.lo)){
    if(add){
      par(mfg=c(2,1))
      xlim.lo = par()$usr[1:2]
    }else if(!is.null(vdir.lo)){
      maxdist = sapply(1:Ns, function(i) max(sapply(1:nrow(vdir.lo), function(j) vdir.lo[j,]%*%x$M[i,,]%*%vdir.lo[j,])))
      # here we must compute the max M projected on vd.lo
      xlim.lo = c(0, max(fk*maxdist))
    }
  }
  if(!is.null(xlim.up)){
    xseq.up = seq(from=xlim.up[1], to=xlim.up[2], length.out=xlength)
  }else{ xseq.up=NULL}
  if(!is.null(xlim.lo)){
    xseq.lo = seq(from=xlim.lo[1], to=xlim.lo[2], length.out=xlength)
  }else{ xseq.lo=NULL}
  
  opar = par()
  opar = par_remove_readonly(opar)
  
  if(closeplot) on.exit(par(opar))
  
  getVdens = function(vdir, xseq){
    if(is.null(vdir)|is.null(xseq)) return(NULL)
    Vdens = sapply(1:nrow(vdir), function(k){
      X = outer(xseq, vdir[k,])
      Y = X[1,,drop=F]*0
      gsi.calcCgram(X,Y,x,FALSE)
    })
    dim(Vdens) = c(Dv,xlength,Dv,nrow(vdir))
    if(!cov){
      # convert to variogram if cov=FALSE
      Y = matrix(rep(0,Dg), ncol=Dg)
      C0 = gsi.calcCgram(Y,Y,x,FALSE)
      dim(C0) = c(Dv,Dv)
      Vdens = sweep(-Vdens, c(1,3), C0, "+") ## this must be corrected when we allow non-symmetric covariances
    }
    return(Vdens)    
  }
  Vdens.up = getVdens(vdir.up, xseq.up)
  Vdens.lo = getVdens(vdir.lo, xseq.lo)
  
  myplot = function(...) matplot(type="l",ylab="", xlab="",xaxt="n", ...)
  if(add) myplot = function(...) matlines(...)
  if(!add){
    par(mfrow=c(Dv+1,Dv+1), mar=c(2,3,0,0), oma=c(1,4,1,1), xpd=NA)
    myplot(c(0,0), c(0,0), pch="", ann=FALSE, bty="n", yaxt="n") 
  }
  for(i in 1:Dv){
      for(j in 1:Dv){
        if((i>=j)&!(is.null(vdir.lo)|is.null(xlim.lo))){
            par(mfg=c(i+1,j,Dv+1,Dv+1))
            ylim = range(Vdens.lo[,,j,])
            if(commonAxis) ylim=range(Vdens.lo[,,j,])
            myplot(xseq.lo, Vdens.lo[i,,j,], ylim=ylim, ...)
            if(i==j & !add){
              axis(side = 3)
              mtext(text=varnames[i], side = 3, line=3)
              mtext(text=varnames[i], side = 4, line=3)
            }
        }
        if((i<=j)&!(is.null(vdir.up)|is.null(xlim.up))){
            par(mfg=c(i,j+1,Dv+1,Dv+1))
            ylim = range(Vdens.up[,,j,])
            if(commonAxis) ylim=range(Vdens.up[,,j,])
            myplot(xseq.up, Vdens.up[i,,j,], ylim=ylim,...)
            if(i==j & !add){
              axis(side = 1)
              mtext(text=varnames[i], side = 1, line=3)
              mtext(text=varnames[i], side = 2, line=3)
            }
        } 
      }}  
  mtext(text="lag distance", side=1, outer = TRUE, line=0)
  mtext(text=c("semivariogram","covariance")[cov+1], side=2, outer = TRUE, line=2)
  invisible(opar)
}




#' Check for anisotropy of a theoretical variogram
#' 
#' Checks for anisotropy of a theoretical variogram or covariance function model

#' @param x variogram or covariance model object
#' @param tol tolerance for 
#' @param ... extra arguments for generic  functionality
#'
#' @return Generic function. Returns of boolean answering the question of the name,
#' or NA if object \code{x} does not contain a known theoretical variogram
#' @export
is.isotropic <- function(x, tol=1e-10, ...){ UseMethod("is.isotropic", x) }

#' @method is.isotropic default
#' @export
is.isotropic.default = function(x, tol=1e-10, ...) NA

#' @method is.isotropic gmCgram
#' @export
is.isotropic.gmCgram = function(x, tol=1e-10, ...){
   all(apply(x$M, 1, function(y){
     ev = eigen(y, only.values=TRUE)[[1]]
     all(abs(ev-ev[1])<tol)
   }))
}

#' @method is.isotropic variogramModel
#' @export
is.isotropic.variogramModel = function(x, tol=1e-10, ...){
  anis = x[,grep("anis", colnames(x))]
  all(apply(anis, 2, function(y) all(abs(y-y[1])<tol) ) )
}

#' @method is.isotropic variogramModelList
#' @export
is.isotropic.variogramModelList = function(x, tol=1e-10, ...) is.isotropic(x[[1]], tol=tol)

#' @method is.isotropic LMCAnisCompo
#' @export
is.isotropic.LMCAnisCompo = function(x, tol=1e-10, ...){
  all(sapply(x["A",], 1, function(y){
    ev = eigen(y$A, only.values=TRUE)[[1]]
    all(abs(ev-ev[1])<tol)
  }))
}









#### empirical variogram ---------------------




#' Variogram method for gmSpatialModel objects
#' 
#' Compute the empirical variogram of the conditioning data contained in a [gmSpatialModel-class] object
#'
#' @param object a gmSpatialModel object containing spatial data.
#' @param methodPars (currently ignored)
#' @param ... further parameters to [gstat::variogram()]
#'
#' @return Currently the function is just a convenience wrapper on
#' the variogram calculation functionalities of package "gstat", 
#' and returns objects of class "\code{gstatVariogram}". Check the
#' help of \code{gstat::variogram} for further information.
#' In the near future, methods will be created, which will depend on
#' the properties of the two arguments provided,  \code{object} and
#' \code{methodPars}.
#' @export
#' @importFrom gstat variogram
variogram_gmSpatialModel <-  function(object, methodPars=NULL, ...){
  if(!is.null(methodPars)) stop("use 'variogram' with named parameters only")
  gstat::variogram(as.gstat(object), ...)
}


# Variogram calculations
# 
# Compute empricial variograms out of a spatial data object
# 
# @param object spatial data container
# @param ... further parameters for variogram calculation
# 
# @return depending on the input data, different kinds of empirical variograms 
# will be produced. See appropriate method descriptions.
# 
# @importFrom sp variogram
# @export
##variogram <- function(object, ...) UseMethod("variogram", object)

# @describeIn variogram
# @method variogram default
# @export
#variogram.default <- function(object, ...){
#  return(variogram_gmSpatialModel(object, ...))
#}

  
#' @export
setGeneric("variogram", function(object,...) standardGeneric("variogram"))  



#' Empirical variogram or covariance function in 2D
#' 
#' compute the empirical variogram or covariance function in a 2D case study 
#' 
#' @param X matrix of Nx2 columns with the geographic coordinates
#' @param Z matrix or data.frame of data with dimension (N,Dv)
#' @param Ff for variogram, matrix of basis functions with nrow(Ff)=N (can be a N-vector of 1s); 
#' for covariance function, a (N,Dv)-matrix or a  Dv-vector giving the mean values
#' @param maxdist maximum lag distance to consider
#' @param lagNr number of lags to consider
#' @param lags if maxdist and lagNr are not specified, either: (a) a matrix of 2 columns giving 
#' minimal and maximal lag distance defining the lag classes to consider, or (b) a vector of lag breaks 
#' @param azimuthNr number of azimuths to consider
#' @param azimuths if azimuthNr is not specified, either: (a) a matrix of 2 columns giving 
#' minimal and maximal azimuth defining the azimuth classes to consider, or (b) a vector of azimuth breaks
#' @param maxbreadth maximal breadth (in lag units) orthogonal to the lag direction
#' @param minpairs minimal number of pairs falling in each class to consider the calculation sufficient; defaults to 10
#' @param cov boolean, is covariance (TRUE) or variogram (FALSE) desired? defaults to variogram
#'
#' @return An empirical variogram for the provided data. NOTE: avoid using directly gsi.* functions! They 
#' represent either internal functions, or preliminary, not fully-tested functions. Use \code{\link{variogram}} instead.
#' @export
#' @family gmEVario functions
#'
#' @examples
#' library(gstat)
#' data("jura", package = "gstat")
#' X = as.matrix(jura.pred[,1:2])
#' Z = as.matrix(jura.pred[,c("Zn","Cd","Pb")])
#' vge = gsi.EVario2D(X,Z)
#' dim(vge)
#' dimnames(vge)
#' class(vge["gamma",1])
#' dim(vge["gamma",1][[1]])
#' vge["npairs",1]
#' vge["lags",1]
gsi.EVario2D = function(X,Z,Ff=rep(1, nrow(X)),
                      maxdist= max(dist(X[sample(nrow(X),min(nrow(X),1000)),]))/2, 
                      lagNr = 15, lags = seq(from=0, to=maxdist, length.out=lagNr+1),
                      azimuthNr=4, azimuths = seq(from=0, to=180, length.out=azimuthNr+1)[1:azimuthNr],
                      maxbreadth=Inf, minpairs=10, cov=FALSE){

  # dimensions
  N = nrow(X)
  Dv = ncol(Z)
  Dg = ncol(X)
  stopifnot(N==nrow(Z))
  if(length(dim(Ff))==0){
    stopifnot(N==length(Ff))
  }else{
    stopifnot(N==nrow(Ff))
  }
  # expand the information given into a set of columns stating conditions
  if(length(dim(lags))==0){
    lags = data.frame(minlag=lags[-length(lags)], maxlag=lags[-1]) 
    if(maxbreadth!=Inf) lags[,"maxbreadth"]=maxbreadth
  }else if(dim(lags)==2){
    lags = data.frame(lags)
    colnames(lags) = c("minlag","maxlag","maxbreadth")[1:ncol(lags)]
  }else stop("lags can be either a vector of lags or a data.frame, see ?gsi.EVario2D")

  if(length(dim(azimuths))==0){
    tol = (azimuths[2]-azimuths[1])/2
    if(is.na(tol)) tol=180
    azimuths = data.frame(minaz=azimuths-tol, maxaz=azimuths+tol) 
  }else if(dim(azimuths)==2){
    azimuths = data.frame(azimuths)
    colnames(azimuths) = c("minaz","maxaz")
  }else stop("azimuths can be either a vector of lags or a data.frame, see ?gsi.EVario2D")
  
  # compute residuals
  if(cov){ 
    kk = 1
    if(dim(Z)==dim(Ff)){
      Z = Z-Ff
    }else if(ncol(Z)==length(c(unlist(Ff)))){
      Z = sweep(Z, 2, Ff, "-")
    }
  }else{
    kk = 1 # 2 # we consider each pair only once
    Z = lm(as.matrix(Z)~as.matrix(Ff)+0)$residuals ## ideally this should be a GLS fit
  }
  # compute pairs
  ij = expand.grid(1:nrow(X), 1:nrow(X))# indices
  op = ifelse(cov, "*", "-")
  ZZ = outer(as.matrix(Z),as.matrix(Z), op) # variables
  ZZ = aperm(ZZ, c(1,3,2,4))
  dim(ZZ) = c(N*N, Dv, Dv)
  XX = X[ij[,1],]-X[ij[,2],] # locations
  XXabs = gmApply(XX, 1, function(x) sqrt(sum(x^2))) 
  XXaz = gmApply(XX, 1, function(x) pi/2-atan2(x[2],x[1])) +2*pi
  XXaz = XXaz %% pi  # residual to 180°
  
  # output

  ## ATTENTION: needs to be changed to return a structure (3,Na)-matrix of objects,
  #     like logratioVariogramAnisotropy
  Nh = nrow(lags)
  Na = nrow(azimuths)
  vg = array(0, dim=c(Nh, Dv, Dv, Na))
  n = array(0, dim=c(Nh, Na))
  azs = azimuths * pi/180
  res = sapply(1:Na, function(i){
    tk_a = (azs[i,1]<=XXaz) & (azs[i,2]>=XXaz)
    zz = ZZ[tk_a,,]
    xxabs = XXabs[tk_a]
    xxaz = XXaz[tk_a]
    tk_h = outer(xxabs, lags[,1],">=") & outer(xxabs, lags[,2],"<=")
    if(ncol(lags)>2){
      tk_b = outer(xxabs * abs(sin((xxaz-(azs[i,2]-azs[i,1])))), lags[,3], "<=")
      tk_h = tk_h & tk_b
    } 
    n[,i] = colSums(tk_h)
    for(j in 1:Nh){
      if(n[j,i]>minpairs){
        vg[j,,,i] = gmApply(zz[tk_h[,j],,], c(2,3),"sum")/(kk*n[j,i])
      }else{
        vg[j,,,i]=NA
      }
    }
    return(list(gamma=vg[,,,i], lags=gsi.lagClass(lags), npairs =n[,i]))
  })
  
  # output
  attr(res, "directions") = gsi.azimuthInterval(azimuths)
  # attr(res, "lags") = gsi.lagClass(lags)
  attr(res, "type") = ifelse(cov, "covariance","semivariogram")
  class(res) = "gmEVario"
  return(res)
}




#' Plot empirical variograms
#' 
#' Flexible plot of an empirical variogram of class gmEVario
#' 
#' @param x object to print, of class gmEVario
#' @param xlim.up range of X values to be used for the diagrams of the upper triangle
#' @param xlim.lo range of X values to be used for the diagrams of the lower triangle
#' @param vdir.up in case of anisotropic variograms, indices of the directions to be plotted 
#' on the upper triangle
#' @param vdir.lo ..., indices of the directions to be plotted on the lower triangle
#' @param varnames variable names to be used
#' @param type  string, controlling whether to plot lines, points, etc (see \code{\link{plot}})
#' @param add boolean, add stuff to an existing diagram?
#' @param commonAxis boolean, should vertical axes be shared by all plots in a row?
#' @param cov boolean, is this a covariance? (if FALSE, it is a variogram)
#' @param closeplot logical, should the plot be left open (FALSE) for further changes, or be frozen (TRUE)? 
#' defaults to TRUE
#' @param ... further parameters to \code{\link{matplot}}
#'
#' @return invisibly, the graphical parameters active before calling the function. 
#' This is useful for freezing the plot if you provided `closeplot=FALSE`.
#' 
#' How to use arguments `vdir.lo` and `vdir.up`? Each empirical variogram \code{x} has been 
#' computed along certain distances, recorded in its attributes and retrievable with command
#' \code{\link{ndirections}}.   
#' @export
#' @family gmEVario functions
#' @method plot gmEVario
#'
#' @examples
#' library(gstat)
#' data("jura", package = "gstat")
#' X = as.matrix(jura.pred[,1:2])
#' Z = as.matrix(jura.pred[,c("Zn","Cd","Pb")])
#' vge = gsi.EVario2D(X,Z)
#' plot(vge)
#' plot(vge, pch=22, lty=1, bg="grey")
plot.gmEVario = function(x, xlim.up=NULL, xlim.lo=NULL, vdir.up= NULL, vdir.lo= NULL, 
                         varnames = dimnames(x$gamma)[[2]], type="o", 
                         add=FALSE, commonAxis=TRUE, cov =attr(x,"type")=="covariance",
                         closeplot=TRUE, ...){
  Dv = dim(x[1,1][[1]])[2]
  if(is.null(varnames)) varnames = paste("v", 1:Dv, sep="")
  if(is.null(vdir.up)&is.null(vdir.lo)) vdir.lo <- 1:ndirections(x)
  if( any(c(vdir.up, vdir.lo)>ndirections(x))){
    stop("indicated directions (vdir.up or vdir.lo) do not exist in x")
  }
  if(is.null(xlim.up)){
    if(add){
      par(mfg=c(1,2))
      xlim.up = par()$usr[1:2]
    }else if(!is.null(vdir.up)){
      maxdist = max(sapply(x["lags",], gsi.midValues.lagClass ) )
      xlim.up = c(0, maxdist)
    }
  }
  if(is.null(xlim.lo)){
    if(add){
      par(mfg=c(2,1))
      xlim.lo = par()$usr[1:2]
    }else if(!is.null(vdir.lo)){
      maxdist = max(sapply(x["lags",], gsi.midValues.lagClass ) )
      xlim.lo = c(0, maxdist)
    }
  }

  opar = par()
  opar = par_remove_readonly(opar)
  
  if(closeplot) on.exit(par(opar))
  
  myplot = function(...) matplot(type=type, ylab="", xlab="",xaxt="n", ...)
  if(add) myplot = function(...) matpoints(type=type, ...)
  if(!add){
    myplot(c(0,0), c(0,0), pch="", ann=FALSE, bty="n", yaxt="n")
    par(mfrow=c(Dv+1,Dv+1), mar=c(2,3,0,0), oma=c(1,4,1,1), xpd=NA)
  }
  for(i in 1:Dv){
    for(j in 1:Dv){
      if((i>=j)&(!is.null(vdir.lo))){
        par(mfg=c(i+1,j,Dv+1,Dv+1))
        ylim = range(sapply(vdir.lo, function(kk) x["gamma",kk][[1]][,i,j]))
        if(commonAxis) ylim=range(sapply(vdir.lo, function(kk) x["gamma",kk][[1]][,,j]))
        myplot(
          sapply(vdir.lo, function(kk) gsi.midValues.lagClass(x["lags",kk][[1]])), 
          sapply(vdir.lo, function(kk) x["gamma",kk][[1]][,i,j]),  ylim=ylim, ...)
        if(i==j){
          axis(side = 3)
          mtext(text=varnames[i], side = 3, line=3)
          mtext(text=varnames[i], side = 4, line=3)
        }
      }
      if((i<=j)&(!is.null(vdir.up))){
        par(mfg=c(i,j+1,Dv+1,Dv+1))
        ylim = range(sapply(vdir.up, function(kk) x["gamma",kk][[1]][,i,j]))
        if(commonAxis) ylim=range(sapply(vdir.up, function(kk) x["gamma",kk][[1]][,,j]))
        myplot(
          sapply(vdir.up, function(kk) gsi.midValues.lagClass(x["lags",kk][[1]])), 
          sapply(vdir.up, function(kk) x["gamma",kk][[1]][,i,j]),  ylim=ylim, ...)
        if(i==j){
          axis(side = 1)
          mtext(text=varnames[i], side = 1, line=3)
          mtext(text=varnames[i], side = 2, line=3)
        }
      } 
    }}  
  mtext(text="lag distance", side=1, outer = TRUE, line=0)
  mtext(text=c("semivariogram","covariance")[cov+1], side=2, outer = TRUE, line=2)
  attr(opar, "vdir.up")=vdir.up
  attr(opar, "vdir.lo")=vdir.lo
  invisible(opar)
}






#' Convert empirical structural function to gmEVario format 
#' 
#' Convert empirical covariance functions or variograms to the format gmEVario 
#' of package gmGeostats
#' 
#' @param vgemp variogram/covariance function to be converted
#' @param ... further parameters
#'
#' @return the empirical covariance function or variogram, recasted to class 
#' \code{gmEVario}. This is a generic function. Methods exist for objects of 
#' class \code{logratioVariogram}\code{logratioVariogramAnisotropy} 
#' (for compositional data) and \code{gstatVariogram}
#' (from package \code{gstat}).  
#' @export
#' @aliases as.gmEVario.gstatVariogram as.gmEVario.logratioVariogram
#'  as.gmEVario.logratioVariogramAnisotropy
#'
#' @family gmEVario functions
as.gmEVario  <- function(vgemp,...){ UseMethod("as.gmEVario",vgemp)}
as.gmEVario.default  <- function(vgemp,...) vgemp




#' @describeIn variogramModelPlot.gmEVario Quick plotting of empirical and theoretical variograms
#' @export
variogramModelPlot <- function(vg, ...)  UseMethod("variogramModelPlot", vg)


#' Quick plotting of empirical and theoretical variograms
#' Quick and dirty plotting of empirical variograms/covariances with or without their models 
#' @param vg empirical variogram or covariance function
#' @param model optional, theoretical variogram or covariance function
#' @param col colors to use for the several directional variograms
#' @param commonAxis boolean, should all plots in a row share the same vertical axis?
#' @param newfig boolean, should a new figure be created? otherwise user should ensure the device space is appropriately managed  
#' @param closeplot logical, should the plot be left open (FALSE) for further changes, or be frozen (TRUE)? 
#' defaults to TRUE
#' @param ... further parameters to underlying plot or matplot functions
#'
#' @return The function is primarily called for producing a plot. However, it 
#' invisibly returns the graphical parameters active before the call 
#' occurred. This is useful for constructing complex diagrams, by adding layers 
#' of info. If you want to "freeze" your plot, embed your call in another
#' call to \code{\link{par}}, e.g. \code{par(variogramModelPlot(...))}; if you
#' want to leave the plot open for further changes give the extra argument `closeplot=FALSE`.
#' @export
#' @family variogramModelPlot 
#' @family gmEVario functions 
#' @family gmCgram functions 
#' @seealso [logratioVariogram()]
#' @method variogramModelPlot gmEVario
#'
#' @examples
#' utils::data("variogramModels")
#' v1 = setCgram(type=vg.Gau, sill=diag(3)+0.5, anisRanges = 5e-1*diag(c(3,0.5)))
#' v2 = setCgram(type=vg.Exp, sill=0.3*diag(3), anisRanges = 5e-2*diag(2))
#' vm = v1+v2
#' plot(vm, closeplot=TRUE)
#' library(gstat)
#' data("jura", package = "gstat")
#' X = as.matrix(jura.pred[,1:2])
#' Z = as.matrix(jura.pred[,c("Zn","Cd","Pb")])
#' vge = gsi.EVario2D(X,Z)
#' variogramModelPlot(vge, vm)
#' 
#' 
variogramModelPlot.gmEVario <- function(vg, model = NULL,   # gstat  or variogramModelList object containing a variogram model fitted to vg
                                        col = rev(rainbow(ndirections(vg))),
                                        commonAxis = FALSE,
                                        newfig = TRUE, closeplot=TRUE, ...){
  opar = plot(vg, commonAxis=commonAxis, add=!newfig, col=col, closeplot=is.null(model), ...)
  opar = par_remove_readonly(opar)
  
  if(closeplot) on.exit(par(opar))
  vdir.lo = attr(opar, "vdir.lo")
  vdir.up = attr(opar, "vdir.up")
  if(is.null(model))
    return(invisible(opar))
  # OTHERWISE: add the curves for the model
  aux = as.directorVector(attr(vg, "directions"))
  if(!is.null(vdir.lo)) vdir.lo = aux[vdir.lo,] 
  if(!is.null(vdir.up)) vdir.up = aux[vdir.up,] 
  opar = plot(as.gmCgram(model), vdir.up= vdir.up, vdir.lo= vdir.lo, add=TRUE, cov =FALSE, ...)
  #f = as.function(as.gmCgram(gg))
  #Dv = dim(vg$gamma)[2]
  #dirs = as.directorVector(attr(vg, "directions"))
  #Dg = ncol(dirs)
  #for(i in 1:Dv){
  #  for(j in 1:Dv){
  #    if((i>=j)&(!is.null(vdir.lo))){
  #      par(mfg=c(i+1,j,Dv+1,Dv+1))
  #      dirs.lo = dirs[vdir.lo,, drop=F]
  #      xlim = par()$usr[1:2]
  #      hd = seq(from=xlim[1], to=xlim[2], length.out = 200)
  #      Y = outer(hd, 1:nrow(dirs.lo), function(h,k){
  #        f(rep(0,Dg), dirs.lo[k,]*h)[i,j]
  #      }) 
  #      matlines(hd, Y, col=col[vdir.lo], ...)
  #    }
  #    if((i<=j)&(!is.null(vdir.up))){
  #      par(mfg=c(i,j+1,Dv+1,Dv+1))
   #     dirs.up = dirs[vdir.up,, drop=F]
  #      xlim = par()$usr[1:2]
  #      hd = seq(from=xlim[1], to=xlim[2], length.out = 200)
  #      Y = outer(hd, 1:nrow(dirs.up), function(h,k){
  #        f(rep(0,Dg), dirs.up[k,]*h)[i,j]
  #      }) 
  #      matlines(hd, Y, col=col[vdir.up], ...)
  #    } 
  #  }}  
  invisible(opar)
} 



#' Number of directions of an empirical variogram
#' 
#' Returns the number of directions at which an empirical variogram was computed
#' 
#' @param x empirical variogram object
#'
#' @return Generic function. It provides the 
#' number of directions at which an empirical variogram was computed
#' @export
#' @family gmEVario functions 
#' @family gmCgram functions 
#' @seealso [logratioVariogram()], [gstat::variogram()]
ndirections <- function(x){ UseMethod("ndirections", x) }

#' @describeIn ndirections generic method
#' @method ndirections default
#' @export
ndirections.default = function(x) nrow(x)
#' @describeIn ndirections method for objects of class "azimuth" (vectors of single angles)
#' @method ndirections azimuth
#' @export
ndirections.azimuth = function(x) length(x)
#' @describeIn ndirections method for objects of class "azimuthInterval" (data.frames of intervals for angles)
#' @method ndirections azimuthInterval
#' @export
ndirections.azimuthInterval = function(x) length(x[[1]])
#' @describeIn ndirections method for empirical logratio variograms with anisotropy 
#' @method ndirections logratioVariogramAnisotropy
#' @export
ndirections.logratioVariogramAnisotropy = function(x) ndirections(attr(x,"directions"))
#' @describeIn ndirections method for empirical logratio variograms without anisotropy
#' @method ndirections logratioVariogram
#' @export
ndirections.logratioVariogram = function(x) 1
#' @describeIn ndirections method for empirical gmGeostats variograms
#' @method ndirections gmEVario
#' @export
ndirections.gmEVario = function(x) ndirections(attr(x,"directions"))
#' @describeIn ndirections method for empirical gstat variograms
#' @method ndirections gstatVariogram
#' @export
ndirections.gstatVariogram = function(x){
  length(unique(paste(x$dir.hor, x$dir.ver)))
} 





#### internal functions -------------

gsi.azimuth = function(x){
   class(x) = c("azimuth","directionClass")
   return(x)
} 

gsi.azimuthInterval = function(x){
  class(x) = c("azimuthInterval","directionClass")
  return(x)
} 



gsi.directorVector = function(x){
  if(length(dim(x))!=2) stop("provided director vectors are not a matrix!")
  class(x) = c("directorVector", "directionClass")
  return(x)
}

print.directionClass = function(x, complete=TRUE, ...){
  cat(paste("with",nrow(as.directorVector(x)),"directions\n"))
  if(complete){
    print(unclass(x), ...)
  }
} 




as.directorVector <- function(x){ UseMethod("as.directorVector",x) }

#' @method as.directorVector default 
as.directorVector.default = function(x,...) x

#' @method as.directorVector azimuth
as.directorVector.azimuth = function(x, D=2){
  res = cbind(cos(pi/2-x), sin(pi/2-x))
  if(D>2){
    res = cbind(res, matrix(0, ncol=D-2, nrow=nrow(res)))
  }
  colnames(res) = paste("v", 1:ncol(res), sep="")
  return(gsi.directorVector(res))
}

#' @method as.directorVector azimuthInterval
as.directorVector.azimuthInterval = function(x, D=2){
  res = (x[[1]]+x[[2]])/2
  return(as.directorVector.azimuth(res))
}


gsi.lagdists = function(x){
  class(x) = c("lagdist","lagClass")
  return(x)
} 

gsi.lagClass = function(x){
  if(length(dim(x))!=2) stop("provided lag classes are not matrix-like!")
  class(x) = c("lagClass")
  return(x)
}

#' @method print lagClass
print.lagClass = function(x,...){
  if(is.list(x)){
    print(as.data.frame(unclass(x)), ...)
  }else{
    print(unclass(x), ...)
  }
}

gsi.midValues <- function(x) UseMethod("gsi.midValues",x)

#' @method gsi.midValues default
gsi.midValues.default = function(x) x

#' @method gsi.midValues lagClass
gsi.midValues.lagClass <- function(x){ (x[[1]]+x[[2]])/2 }
  
#' @method gsi.midValues azimuthIntervals
gsi.midValues.azimuthInterval <- function(x){ (x[[1]]+x[[2]])/2 }





## theoretical structural functions
# S3 -> S4 classes  
# cat("creating variogram model classes\n")
#' @include gstatCompatibility.R
#' @include compositionsCompatibility.R
setOldClass("gmCgram")
setOldClass("LMCAnisCompo")
setOldClass("variogramModelList")
setOldClass("variogramModel")

# abstract classes
#' @title Structural function model specification
#' @description Abstract class, containing any specification of a variogram (or covariance) model
#' @export
setClassUnion(name="ModelStructuralFunctionSpecification", 
              members=c("NULL","gmCgram", "LMCAnisCompo", "variogramModelList", "variogramModel"))



# 
# #### container class --------------
# # An S4 class to represent a Gaussian random field specification
# #
# # @slot structure ModelStructuralFunctionSpecification. Variogram or 
# # (generalised) covariance function specification, typically an object
# # obtained from a call to functions such as \code{\link{setCgram}},
# # \code{\link{LMCAnisCompo}} or \code{gstat::vgm}.
# # @slot formula formula specifying the structure
# # of dependence of the mean of the random field w.r.to spatial coordinates 
# # and/or covariables; typically it will have no left-hand-side term; 
# # @slot beta numeric, a vector with as many coefficients as terms the formula
# # above requires for a full specification of the trend; if unknown, these can
# # be NAs, as many as needed. 
# #
# # @return A object with the slots populated as given
# # @export
# # @seealso [gmSpatialModel-class], and the `make.gm*` functions referenced there
# setClass("gmGaussianModel", 
#          slots = list(structure = "ModelStructuralFunctionSpecification",
#                       formula="formula",
#                       beta = "structure")
# )
# 
# #setMethod("initialize", signature="gmGaussianModel",
# #          def=function(.Object, structure, formula, beta){
# #            .Object@formula = formula
# #            .Object@beta = beta
# #            if(!is.null(structure)) .Object@structure = structure
# #            return(.Object)
# #          }
# #)


#### container class --------------
# An S4 class to represent a Gaussian random field specification
#
# @slot structure ModelStructuralFunctionSpecification. Variogram or 
# (generalised) covariance function specification, typically an object
# obtained from a call to functions such as \code{\link{setCgram}},
# \code{\link{LMCAnisCompo}} or \code{gstat::vgm}.
# @slot formula formula specifying the structure
# of dependence of the mean of the random field w.r.to spatial coordinates 
# and/or covariables; typically it will have no left-hand-side term; 
# @slot beta numeric, a vector with as many coefficients as terms the formula
# above requires for a full specification of the trend; if unknown, these can
# be NAs, as many as needed. 
#
# @return A object with the slots populated as given
# @export
# @seealso [gmSpatialModel-class], and the `make.gm*` functions referenced there
setClass("gmGaussianModel", 
         slots = list(structure = "ModelStructuralFunctionSpecification",
                      formula="formula",
                      beta = "numeric")
)


## empirical structural functions
# S3 -> S4 classes
# cat("creating empirical variogram classes\n")
setOldClass("gmEVario")
setOldClass("logratioVariogram")
setOldClass("logratioVariogramAnisotropy")
setOldClass("gstatVariogram")


# abstract classes
#' @title Empirical structural function specification
#' @description Abstract class, containing any specification of an empirical variogram (or covariance function, or variations)
#' @export
setClassUnion(name="EmpiricalStructuralFunctionSpecification", members=c("NULL","gmEVario", "logratioVariogram", "logratioVariogramAnisotropy", "gstatVariogram"))

