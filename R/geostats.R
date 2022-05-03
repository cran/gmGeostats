# library("compositions")

# ATTENTION: when calling this script, it might be required to 
#     define compileme=TRUE or FALSE, as well as to copy/link
#     this directory as a subdirectory of the working directory
# 
# if( Sys.info()["sysname"]=="Linux"){
#   if(!exists("compileme")) compileme =FALSE
#   if(compileme) system("R CMD SHLIB -fopenmp src/CGeostats.c")
# 
#   source("R/spSupport.R")
#   source("R/variograms.R")
#   source("R/Anamorphosis.R")
# 
#   if( is.loaded("CMVTurningBands" )) {
#     dyn.unload("src/CGeostats.so")
#   }
#   dyn.load("src/CGeostats.so")
# }else{
#   source("R/spSupport.R")
#   source("R/variograms.R")
#   source("R/Anamorphosis.R")
#  # if( is.loaded("CMVTurningBands" )) {
#  #   dyn.unload("src/gmGeostats.dll")
#  # }
# #  dyn.load("src/gmGeostats.dll")
#   }


## 
#' Compute covariance matrix oout of locations 
#'
#' internal function to compute the variogram model for all combinations of two 
#     sets of coordinates.
#'
#' @param X matrix, coordinates of the first set of locations 
#' @param Y matrix, coordinates of the first set of locations (often Y=X, and ijEqual=TR)
#' @param vgram covariogram model, of format "gmCgram"
#' @param ijEqual logical, are X==Y? if you put TRUE, the function will only compute half of
#' the points, so be sure that this is right!
#'
#' @return a matrix of elements of the covariance, in NxM blocks of DxD elements
#' (being N and M resp. the nr. of rows of X and Y, and D the nr of variables) 
#' @useDynLib gmGeostats CcalcCgram
gsi.calcCgram <- function(X,Y,vgram,ijEqual=FALSE) {
  d <- nrow(vgram$nugget)
  m <- ncol(X)
  stopifnot(m==ncol(Y))
  nX<-nrow(X)
  nY<-nrow(Y)
  k <- length(vgram$type)
  dimC = c(d,nX,d,nY)
  erg<- .C("CcalcCgram",
           dimX=checkInt(dim(X),2),
           ldX=checkInt(dim(X)[1],1),
           X=checkDouble(X),
           dimY=checkInt(dim(Y),2),
           ldY=checkInt(dim(Y)[1],1),
           Y=checkDouble(Y),
           dimC=checkInt(dimC),
           C=numeric(prod(dimC)),
           Nugget=checkDouble(vgram$nugget,d*d),
           nCgrams=checkInt(k,1),
           typeCgrams=checkInt(vgram$type,k),
           A = checkDouble(t(apply(vgram$M,1,gsi.powM,alpha=-1/2)),k*m*m),
           Sill=checkDouble(vgram$sill,k*d*d),
           moreCgramData=checkDouble(vgram$data,k),
           ijEqual=checkLogical(ijEqual,1)#,     PACKAGE = "gmGeostats"
  )
  structure(erg$C,dim=c(d*nX,d*nY))
}

gsiGetUnitVec <- function(dimX,ip) {
    erg = .C("getUnitvecR",dimX=checkInt(dimX,1),ip=checkInt(ip,1),unitvec=numeric(dimX))
    unitvec = erg$unitvec
    unitvec
    }



## 
#' Internal function, unconditional turning bands realisations
#' 
#' Internal function to compute unconditional turning bands simulations
#'
#' @param X matrix of coordinates of locations where realisations are desired
#' @param vgram covariogram model, of format "gmCgram"
#' @param nBands number of bands to use
#' @param nsim number of realisations to return
#' 
#' @return an array with (npoint, nvar, nsim)-elements, being npoint=nrow(X)
#' and nvar = nr of variables in vgram
#' @export
#' @useDynLib gmGeostats CMVTurningBands
gsi.TurningBands <- function(X,vgram,nBands,nsim=NULL) {
  # checks and preps
  stopifnot(length(as.integer(nBands))==1)
  d <- dim(vgram$nugget)[1]
  m <- ncol(X)
  k<-length(vgram$type)
  stopifnot( all(dim(vgram$nugget)==c(d,d)))
  stopifnot( all(dim(vgram$sill)==c(k,d,d)))
  stopifnot( all(dim(vgram$M)==c(k,m,m)))
  # extend to 3D if 2D
  m0 = 3
  if(m<3){
    m0 = m
    X = cbind(X, 0, 0)
    paveId = function(i){
      M = diag(3)
      M[1:m, 1:m] = vgram$M[i,,]
      return(M)
    }
    M = sapply(1:k, paveId)
    dim(M) = c(3,3,k)
    M = aperm(M, c(1,2,3))
    vgram$M = M
    m = 3
  }
  # run!
  if( is.null(nsim) ) {
    erg<-.C("CMVTurningBands",
            dimX=checkInt(dim(X),2),
            X=checkDouble(X,prod(dim(X))),
            dimZ=checkInt(c(d,nrow(X),1),3),
            Z = numeric(nrow(X)*d),
            nBands=checkInt(nBands,1),
            sqrtNugget=checkDouble(gsi.powM(vgram$nugget,1/2),d*d),
            nCgrams = checkInt(length(vgram$type),1),
            typeCgram = checkInt(vgram$type,k),
            A = checkDouble(t(apply(vgram$M,1,gsi.powM,alpha=-1/2)),k*m*m),
            sqrtSill=checkDouble(t(apply(vgram$sill,1,gsi.powM,alpha=1/2)),k*d*d),
            moreCgramData=checkDouble(vgram$data,k)#,     PACKAGE = "gmGeostats"
    )
      erg = cbind(X[, 1:m0],Z=t(structure(erg$Z,dim=c(d,nrow(X)))))
      return(erg) 
  } else {
    nsim <- checkInt(nsim,1)
    stopifnot(nsim>0)
    erg<-.C("CMVTurningBands",
            dimX=checkInt(dim(X),2),
            X=checkDouble(X,prod(dim(X))),
            dimZ=checkInt(c(d,nrow(X),nsim),3),
            Z = numeric(nrow(X)*d*nsim),
            nBands=checkInt(nBands,1),
            sqrtNugget=checkDouble(gsi.powM(vgram$nugget,1/2),d*d),
            nCgrams = checkInt(length(vgram$type),1),
            typeCgram = checkInt(vgram$type,k),
            A = checkDouble(t(apply(vgram$M,1,gsi.powM,alpha=-1/2)),k*m*m),
            sqrtSill=checkDouble(t(apply(vgram$sill,1,gsi.powM,alpha=1/2)),k*d*d),
            moreCgramData=checkDouble(vgram$data,k)#,     PACKAGE = "gmGeostats"
    )
    erg = aperm(structure(erg$Z,dim=c(d,nrow(X),nsim)),c(2,1,3))
    return(erg) 
  }
}





#' Internal function, conditional turning bands realisations
#' 
#' Internal function to compute conditional turning bands simulations
#'
#' @param Xout matrix of coordinates of locations where realisations are desired
#' @param Xin matrix of coordinates of locations of conditioning points
#' @param Zin matrix of variables at the conditioning locations
#' @param vgram covariogram model, of format "gmCgram"
#' @param nbands number of bands to use
#' @param tol tolerance for the inversion of the cokriging matrix (in case of near-singularity)
#' @param nsim number of realisations to return
#'
#' @return an array with (npoint, nvar, nsim)-elements, being npoint=nrow(X)
#' and nvar = nr of variables in vgram
#' @useDynLib gmGeostats CCondSim
gsi.CondTurningBands <- function(Xout, Xin, Zin, vgram, 
                                 nbands=400, tol=1E-15, nsim=NULL) {
  CII <- gsi.calcCgram(Xin,Xin,vgram,ijEqual=TRUE)
  Cinv <- gsiInv(CII,tol=tol)
  d <- dim(vgram$nugget)[1]
  m <- ncol(Xin)
  k<-length(vgram$type)
  nin<-nrow(Xin)
  Zin<-t(Zin)
  X <- rbind(Xin,Xout)
  dimZ <-c(nrow(Zin),nrow(X))
  dimZin <-c(nrow(Zin),nrow(Xin))
  if( is.null(nsim) ) {
    normal <- TRUE
    nsim<-1
  } else {
    normal<-FALSE
  }
  nsim <- checkInt(nsim,1)
  # extend to 3D if 2D
  m0 = 3
  if(m<3){
    m0 = m
    X = cbind(X, 0, 0)
    paveId = function(i){
      M = diag(3)
      M[1:m, 1:m] = vgram$M[i,,]
      return(M)
    }
    M = sapply(1:k, paveId)
    dim(M) = c(3,3,k)
    M = aperm(M, c(1,2,3))
    vgram$M = M
    m = 3
  }
  erg <- .C("CCondSim",
            dimZin=checkInt(dimZin,2),
            Zin   =checkDouble(Zin),
            Cinv  =checkDouble(Cinv,length(Zin)^2),
            dimX  =checkInt(dim(X),2),
            X=checkDouble(X),
            dimZ=checkInt(c(dimZ,nsim),3),
            Z=numeric(prod(dimZ)*nsim),
            nBands=checkInt(nbands,1),
            sqrtNugget=checkDouble(gsi.powM(vgram$nugget,1/2),d*d),
            Nugget=checkDouble(vgram$nugget,d*d),
            nCgrams = checkInt(length(vgram$type),1),
            typeCgram = checkInt(vgram$type,k),
            A = checkDouble(t(apply(vgram$M,1,gsi.powM,alpha=-1/2)),k*m*m),
            sqrtSill=checkDouble(t(apply(vgram$sill,1,gsi.powM,alpha=1/2)),k*d*d),
            Sill=checkDouble(vgram$sill,k*d*d),
            moreCgramData=checkDouble(vgram$data,k),
            cbuf=numeric(d*d*nin),
            dbuf=numeric(d*nin*nsim)#,     PACKAGE = "gmGeostats"
  )
  if( normal ) {
    cbind(Xout,t(structure(erg$Z,dim=dimZ))[-(1:nrow(Xin)),])
  } else {
    aperm(structure(erg$Z,dim=c(dimZ,nsim))[,-(1:nrow(Xin)),,drop=FALSE],c(2,1,3))
  }
}



#### cokriging --------------


### tests --
# if(!exists("do.test")) do.test=FALSE
# if(do.test){
  # library("compositions")
  # library("sp")
  # library("magrittr")
  # library("gstat")
  # data("jura", package = "gstat")
  # dt = jura.pred %>% dplyr::select(Cd:Cr)
  # X = jura.pred[,1:2]
  # spc = SpatialPointsDataFrame(SpatialPoints(X),acomp(dt))
  # a2 = make.gmCompositionalGaussianSpatialModel(spc)
  # #attach(gsi)
  # a3gs = as.gstat(a2)
  # a3vg = variogram(a3gs)
  # a3gs = fit_lmc(a3vg, a3gs, vgm("Sph", psill=1, nugget=0, range=1))
  # a4 = make.gmCompositionalGaussianSpatialModel(spc, V="alr", model=a3gs$model)
  # a4gs = as.gstat(a4)
  # p1 = predict(a4gs, newdata = jura.grid)
  # p2 = predict(a4gs, newdata = jura.grid[,1:2])
  # p3 = predict(a4gs, newdata = SpatialPixels(SpatialPoints(jura.grid[,1:2])))
#   
# }


## 
# Xin, Zin input (coordinates, variables)
# Xout output coordinates
# vgram variogram object
# Fin, Fout: trend basis functions evaluated on each Xin resp. Xout (UK, KT); 
#            default to vectors of 1 (OK); if set to NULL, SK will be used
# krgVar: should the cokriging variance be reported? FALSE (no), TRUE (yes, as a list of matrices
#            containing one entry for each point), NA (yes, as a huge matrix with all points)
# tol: tolerance to determination of quasi-zero eigenvalues in the inversion process
#' Cokriging of all sorts, internal function
#' 
#' internal function to compute cokriging (simple, ordinary, universal, with trend)
#'
#' @param Xout (M,g)-matrix with coordinates of desired interpolation locations
#' @param Xin (N,g)-matrix with coordinates of conditioning locations
#' @param Zin (N,D)-matrix with condtioning observations
#' @param vgram D-dimensional variogram model of class "gmCgram"
#' @param Fin either NULL or a (N,p)-matrix with the base functions of the trend evaluated at 
#' the conditioning locations (defaults to a vector of $N$ 1's)
#' @param Fout either NULL or a (m,p)-matrix with the base functions of the trend evaluated at 
#' the conditioning locations (defaults to a vector of $M$ 1's)
#' @param krigVar logical or NA, should kriging variances be returned? FALSE=no, 
#' TRUE=give point-wise cokriging covariance; NA=return the whole covariance matrix of all points
#' @param tol tolerance for the inversion of the cokriging matrix (in case of near-singularity)
#'
#' @return A (M,D)-matrix of predictions, eventually with an attribute "krigVar"
#' containing the output defined by argument `krigVar`
gsi.Cokriging <- function(Xout,Xin,Zin,vgram,Fin=rep(1,nrow(Xin)),
                          Fout=rep(1,nrow(Xout)),krigVar=FALSE,tol=1E-15){
  stopifnot(is.null(Fin)==is.null(Fout)) # do not continue if Fin and Fout are not compatible
  CII <- gsi.calcCgram(Xin,Xin,vgram,ijEqual=TRUE)
  D = ncol(Zin)
  Id <- diag(D)
  if(!is.null(Fin)){
    FII <- kronecker(Fin, Id)
    CII <- rbind(cbind(CII, FII), cbind(t(FII), 0*Id))
  }
  Cinv <- gsiInv(CII,tol=tol)
  CI <-  gsi.calcCgram(Xin,Xout,vgram,ijEqual=FALSE)
  if(!is.null(Fin)){
    FI <- kronecker(Fout, Id)
    CI <- rbind(CI, t(FI))
  }  
  L <- Cinv %*% CI
  Zout <- t(L)[,1:length(Zin)] %*% c(t(Zin))
  dim(Zout) = c(ncol(Zin), nrow(Xout))
  Zout = t(Zout)
  colnames(Zout) = colnames(Zin)
  if(is.na(krigVar)){
    KV = t(CI) %*% L
  }else{
    if(!krigVar) return(Zout)
    KV = lapply(1:nrow(Xout), function(i){
      idx = (i-1)*D+(1:D)
      t(CI[,idx]) %*% L[,idx]
    })
  }
  attr(Zout,"krigVar") = KV
  return(Zout)
}


