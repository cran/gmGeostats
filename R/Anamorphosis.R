# X: data matrix  ${\bf x}_i, i=1, 2, \ldots, N$
# t: time coordinate
#source("R/spSupport.R")

# if( Sys.info()["sysname"]=="Linux"){
#   if(!exists("compileme")) compileme =FALSE
#   if(compileme) system("R CMD SHLIB -fopenmp ../src/morph.c")
#   
#   if( is.loaded("anaForwardC" )) {
#     dyn.unload("src/morph.so")
#   }
#   dyn.load("src/morph.so")
# }else{
#   if( is.loaded("anaForwardC" )) {
#     dyn.unload("src/gmGeostats.dll")
#   }
#   dyn.load("src/gmGeostats.dll")
# }




################### auxiliary functions ------------------------------
# anaXt: position of each datum after a time "t": ${\bf x}_i(t)$
anaXt <- function(t,X) {
  (1-t)*X
}
# this is a scaling of each row towards the 0-vector: 


# anaSigmat: kernel bandwidth at time "t": $\sigma(t)$
anaSigmat <- function(t,sigma0=0.1,sigma1=1) {
  (sigma1-sigma0)*t+sigma0
}
# this goes from an initial value to the standard normal distribution


# x: an arbitrary point
# X: one of the "pivots" (data)
# t: time instant
# calculate the normal score of the point w.r.t. the sigmat and Xt
anaScore <- function(x,X,t,...){
  (x-anaXt(t,X))/anaSigmat(t,...)
}
  
exp2 = function(x) exp(-0.5*x^2)

# weights: matrix of weights (rows as many as rows of x, columns as many as rows of Y); by rows sum to one
anaW <- function(x,t,Y,sigma0=sigma0,sigma1=sigma1, funkernel = exp2) {
#   sigma = anaSigmat(t,sigma0=sigma0,sigma1=sigma1)
#   W<-matrix(exp(-0.5*c((x[rep(1:nrow(x),nrow(Y)),]-anaXt(t,Y)[rep(1:nrow(Y),each=nrow(x)),])/sigma)^2 %*% rep(1,ncol(Y))),nrow=nrow(x))
  W<-matrix(funkernel(c(((x[rep(1:nrow(x),nrow(Y)),]-anaXt(t,Y)[rep(1:nrow(Y),each=nrow(x)),])/anaSigmat(t,sigma0=sigma0,sigma1=sigma1))^2%*%rep(1,ncol(Y)))),nrow=nrow(x))
  e<- W/c(W%*%rep(1,ncol(W)))
  if(any(e)){
    warning("sum((e%*%rep(1,ncol(W))-1)^2)<0.0001 not satisfied")
  }
  e
}
# x: points that are being deformed
# Y: data points 


# average speed of deformation of the space at time t
anaV <- function(x,t,Y,sigma0,sigma1=1+sigma0) {
  W <- anaW(x,t,Y,sigma0=sigma0,sigma1=sigma1)
  i <- rep(1:nrow(x),each=nrow(Y))
  j <- rep(1:nrow(Y),nrow(x))
  Vind <- matrix(c(t(W))*(-Y[j,]+((sigma1-sigma0)/anaSigmat(t,sigma0=sigma0,sigma1=sigma1))*(x[i,]-anaXt(t,Y)[j,])),nrow=nrow(Y))
  Vall <- matrix(rep(1,nrow(Y))%*%Vind,nrow=nrow(x))
  Vall
}




################### transformation functions ------------------------------


# # forward 
# anaForwardOld <- function(x,Y,sigma0,sigma1=1+sigma0,steps=30,plt=FALSE,sphere=TRUE) {
#   if( sphere ) st = sphTrans(Y) else st=function(x,...) x
#   x<-st(x)
#   Y<-st(Y)
#   h=1/steps
#   if( plt )
#     plot(x[,1:2],pch=".")
#   for(i in 1:steps) {
#     # Verfahren 2.Ordnung
#     v1 = anaV(x,(i-1)*h,Y,sigma0=sigma0,sigma1=sigma1)
#     v2 = anaV(x+h*v1,i*h,Y,sigma0=sigma0,sigma1=sigma1)
#     x  = x+0.5*h*(v1+v2)
#     #x = x+ h*anaV(x,i*h,Y,sigma0=sigma0,sigma1=sigma1)
#     if(plt)
#       points(x[,1:2],pch=".")
#   }
#   x
# }

#' Forward gaussian anamorphosis
#' forward transformation to multivariate gaussian scores
#'
#' @param x points to be transformed (a matrix)
#' @param Y node points defining the transformation (another matrix, same nr. of columns as `x`)
#' @param sigma0 starting spread of the kernels
#' @param sigma1 final spread of the kernels
#' @param steps number of steps to linearize the transform (default 30 is good)
#' @param plt boolean, do you want to get a plot of the transformation?
#' @param sphere boolean, should the data be pre-Y-spherified first? defaults to true
#' @param weights vector of weights for all computations, length must be equal 
#' to number of rows of \code{x}
#'
#' @return a matrix with the gaussian scores; same dimensions of \code{x}
#' @export
#' @seealso [ana()] for defining a function that carries over the transformation  
#' (by means of a closure), [anaBackward()] for the explicit back-transformation, 
#' [sphTrans()] for defining a function that carries over the spherification of the data
#' @useDynLib gmGeostats anaForwardC
#' @author K. Gerald van den Boogaart, Raimon Tolosana-Delgado
#' @examples
#' data("jura", package="gstat")
#' Y = jura.pred[,c(10,12,13)]
#' plot(compositions::acomp(Y))
#' Ylr = compositions::alr(Y)
#' plot(Ylr)
#' z = anaForward(x=Ylr, Y=Ylr, sigma0=0.1)
#' plot(z, asp=1)
#' shapiro.test(z[,1])
#' shapiro.test(z[,2])
anaForward <- function(x,Y,sigma0,sigma1=1+sigma0,steps=30,plt=FALSE,sphere=TRUE, weights=NULL) {
  if(is.null(weights)){
    weights = rep(1,nrow(Y))
  }
  if(length(weights)!=nrow(Y)){
    stop("weights provided not compatible with nrow(Y)")
  }
  # if( sphere ) st = sphTrans(Y) else 
  st=function(x,...) x
  if(is.logical(sphere)){
    if( sphere ) st = sphTrans(Y, weights=weights) else st=function(x,...) x
  }else if(is.function(sphere)){
    st = sphere
  }
  x<-t(st(x))
  Y<-t(st(Y))
  #h=1/steps
  #if( plt )
  #    plot(x[,1:2],pch=".")
  erg<-.C("anaForwardC",
     dimX=checkInt(dim(x),2),
     x=checkDouble(x),
     dimY=checkInt(dim(Y),2),
     y=checkDouble(Y),
     wY=checkDouble(weights,ncol(Y)),
     steps=checkInt(steps,1),
     sigma0=checkDouble(sigma0,1),
     sigma1=checkDouble(sigma1,1),
     PACKAGE = "gmGeostats"
     )
  x<- t(structure(erg$x, dim=dim(x)))
  colnames(x) = paste("flow", 1:ncol(x), sep="")
  return(x)
}


# anaBackwardOld <- function(x,Y,sigma0,sigma1=1+sigma0,steps=30,plt=FALSE,sphere=TRUE) {
#   if( sphere ) st = sphTrans(Y) else st=function(x,...) x
#   Y<-st(Y)
#   h=1/steps
#   if( plt )
#     plot(x[,1:2],pch=".")
#   for(i in 1:steps) {
#     # Verfahren 2.Ordnung
#     v1 = anaV(x,1-(i-1)*h,Y,sigma0=sigma0,sigma1=sigma1)
#     v2 = anaV(x-h*v1,1-i*h,Y,sigma0=sigma0,sigma1=sigma1)
#     x  = x-0.5*h*(v1+v2)
# #    x = x- h*anaV(x,1-i*h,Y,sigma0=sigma0,sigma1=sigma1)
#     if(plt)
#       points(x[,1:2],pch=".")
#   }
#   st(x,inv=TRUE)
# }


#' Backward gaussian anamorphosis
#' backward transformation to multivariate gaussian scores
#'
#' @param x matrix of gaussian scores to be back-transformed
#' @param Y node points defining the transformation (a matrix, same nr of columns)
#' @param sigma0 starting spread of the kernels in the forward transform
#' @param sigma1 final spread of the kernels in the forward transform
#' @param steps number of steps to linearize the transform (default 30 is good)
#' @param plt boolean, do you want to get a plot of the transformation?
#' @param sphere boolean, should the data be taken as pre-Y-spherified? defaults to true
#' @param weights vector of weights for all computations, length must be equal 
#' to number of rows of \code{x}
#'
#' @return a matrix with the scores back-transformed to the same scale as \code{Y}; same dimensions of \code{x}
#' @seealso [ana()] for defining a function that carries over the transformation  
#' (by means of a closure), [anaBackward()] for the explicit back-transformation, 
#' [sphTrans()] for defining a function that carries over the spherification of the data
#' @export
#' @useDynLib gmGeostats anaBackwardC
#' @author K. Gerald van den Boogaart, Raimon Tolosana-Delgado
#' @examples
#' data("jura", package="gstat")
#' Y = jura.pred[,c(10,12,13)]
#' plot(compositions::acomp(Y))
#' Ylr = compositions::alr(Y)
#' Xns = matrix(rnorm(500), ncol=2)
#' plot(Ylr)
#' points(Xns, col=2, pch=4)
#' Xlr = anaBackward(x=Xns, Y=Ylr, sigma0=0.1)
#' qqplot(Xlr[,1], Ylr[,1])
#' qqplot(Xlr[,2], Ylr[,2])
#' qqplot(Xlr[,1]+Xlr[,2], Ylr[,1]+Ylr[,2])
anaBackward <- function(x,Y,sigma0,sigma1=1+sigma0,steps=30,plt=FALSE,sphere=TRUE, weights=NULL) {
  if(is.null(weights)){
    weights = rep(1,nrow(Y))
  }
  if(length(weights)!=nrow(Y)){
    stop("weights provided not compatible with nrow(Y)")
  }
  #if( sphere ) st = sphTrans(Y) else 
  st=function(x,...) x
  if(is.logical(sphere)){
    if( sphere ) st = sphTrans(Y, weights=weights) else st=function(x,...) x
  }else if(is.function(sphere)){
    st = sphere
  }
  Y<-t(st(Y))
  x <- t(x)
  erg<-.C("anaBackwardC",
     dimX=checkInt(dim(x),2),
     x=checkDouble(x),
     dimY=checkInt(dim(Y),2),
     y=checkDouble(Y),
     wY=checkDouble(weights,ncol(Y)),
     steps=checkInt(steps,1),
     sigma0=checkDouble(sigma0,1),
     sigma1=checkDouble(sigma1,1),
     PACKAGE = "gmGeostats"
     )
  x<- t(structure(erg$x, dim=dim(x)))
  st(x,inv=TRUE)
}


# pre-spherification of the data (through SVD)

#' Spherifying transform
#' Compute a transformation that spherifies a certain data set
#'
#' @param Y data set defining the spherifization
#' @param ... extra arguments for generic functionality
#'
#' @return a function with arguments \code{(x, inv=FALSE)}, where \code{x} will be the 
#' data to apply the transformation to, and \code{inv=FALSE} will indicate if the direct 
#' or the inverse transformation is desired. 
#' This function applied to the same data returns a translated, rotated and scaled, so that
#' the new scores are centered, have variance 1, and no correlation.
#' @export
#' @importFrom utils methods
#' @importFrom stats cov.wt
#' 
#' @aliases sphTrans.default
#' @seealso ana, anaBackward, sphTrans
#' @author K. Gerald van den Boogaart, Raimon Tolosana-Delgado
#' 
#' @examples
#' library(compositions)
#' data("jura", package="gstat")
#' Y = acomp(jura.pred[,c(10,12,13)])
#' par(mfrow=c(1,1))
#' plot(Y)
#' sph = sphTrans(Y)
#' class(sph)
#' z = sph(Y)
#' plot(z)
#' cor(cbind(z, ilr(Y)))
#' colMeans(cbind(z, ilr(Y)))
sphTrans <- function(Y,...) UseMethod("sphTrans",Y)

#' @describeIn sphTrans Spherifying transform
#' @method sphTrans default
#' @export
#' @param weights weights to incorporate in the compuations, length=nrow(Y)
#' @param p dimensions to be considered structural (useful for filtering noise)
sphTrans.default <- function(Y, weights=NULL, p=1:ncol(Y), ...) {
  cY = class(Y)[1]
  existsIdt = sub("idt.", "", as.character(methods("idt")))
  if(is.null(weights)){
    if("data.frame" %in% cY){
      b <- gmApply(Y,2,mean)
    }else if(length(intersect(cY, existsIdt))>0){
      b <- mean(idt(Y))
    }else{
      b <- gmApply(Y,2,mean)
    }
    vY = var(idt(Y))
  }else if(length(weights)==nrow(Y)){
      aux = cov.wt(idt(Y),wt=weights)
      b = aux$center
      vY = aux$cov
  }else{
      stop("meaningless weights provided, wrong length?")
  }
  p = p[(p<=ncol(vY))]
  SVD <- svd(vY)
  A <- with(SVD, v[,p] %*% diag(1/sqrt(d[p])) %*% t(u[,p]))
  Ai <- with(SVD, u[,p] %*% diag(sqrt(d[p])) %*% t(v[,p]))
  f<-function(x,inv=FALSE) {
    if(isTRUE(inv)){
      if(class(x)[1]=="matrix"){
        out = t((Ai%*%t(x))+unclass(b))
      }else{
        out = Ai%*%x+b
      }
      return(idtInv(out, orig=Y))
    }else{
      if(class(x)[1]=="matrix"){
        out = t(A%*%(t(x)-unclass(b)))
      }else{
        out = unclass(A %*% (idt(x)-b))
        attr(out, "orig") = NULL
        attr(out, "V") = NULL
      }
      return(out)
    }
  }
  return(f)
}


#' Flow anamorphosis transform
#' Compute a transformation that gaussianizes a certain data set
#'
#' @param Y data set defining the gaussianization
#' @param sigma0 starting spread of the kernels
#' @param sigma1 final spread of the kernels
#' @param steps number of steps to linearize the transform (default 30 is good)
#' @param sphere boolean, should the transform include a spherifization step, 
#' making \code{Y} spherical?
#' @param weights weights to incorporate in the compuations, length=nrow(Y)
#'
#' @return a function with arguments \code{(x, inv=FALSE)}, where \code{x} will be the 
#' data to apply the transformation to, and \code{inv=FALSE} will indicate if the direct 
#' or the inverse transformation is desired 
#' @export
#' @seealso anaForward, anaBackward, sphTrans
#' @author K. Gerald van den Boogaart
#' @examples
#' library(compositions)
#' data("jura", package="gstat")
#' Y = acomp(jura.pred[,c(10,12,13)])
#' plot(Y)
#' anafun = ana(Y)
#' class(anafun)
#' z = anafun(Y)
#' plot(z)
#' y = anafun(z, inv=TRUE)
#' plot(data.frame(orig=Y,recalc=y))
ana <-function(Y,sigma0=0.1,sigma1=1,steps=30,sphere=TRUE, weights=NULL) {
  f<-function(x,inv=FALSE,...){
    if( isTRUE(inv) ) {
      anaBackward(x,Y,sigma0=sigma0,sigma1=sigma1,sphere=sphere,weights=weights,...)
    } else {
      anaForward(x,Y,sigma0=sigma0,sigma1=sigma1,sphere=sphere,weights=weights,...)
    }
  }
  f
}




