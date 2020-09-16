### a test of lack of spatial correlation
#' Test for lack of spatial correlation
#' 
#' Permutation test for checking lack of spatial correlation.
#'
#' @param Z matrix (or equivalent) of scaled observations
#' @param ... extra arguments for generic functionality
#'
#' @return Produces a test of lack of spatial correlation by means of permutations. The test
#' statistic is based on the smallest eigenvalue of the generalised eigenvalues of the matrices
#' of covariance for short range and for long range.
#' @export
#' @importFrom boot boot
#' @examples
#' \donttest{ 
#' data("jura", package="gstat")
#' X = jura.pred[, 1:2]
#' Z = data.frame(compositions::ilr(jura.pred[,-(1:6)]))
#' noSpatCorr.test(Z=Z, X=X)
#' # now destroy the spatial structure reshuffling the coordinates:
#' ip = sample(nrow(X))
#' noSpatCorr.test(Z=Z, X=X[ip,]) 
#' }
noSpatCorr.test <- function(Z, ...) UseMethod("noSpatCorr.test", Z)


#' @describeIn noSpatCorr.test  Test for lack of spatial correlation
#' @method noSpatCorr.test data.frame
#' @param X matrix (or equivalent) of sample location coordinates
#' @export
noSpatCorr.test.data.frame <- function(Z, X, ...){
  noSpatCorr.test(as.matrix(Z), X = X, ...)
}

#' @describeIn noSpatCorr.test  Test for lack of spatial correlation, works only for 
#' Spatial objects with a "data" slot
#' @method noSpatCorr.test default
#' @export
noSpatCorr.test.default <- function(Z, ...){
  if(is(Z, "Spatial")){
    coords = sp::coordinates(Z)
    if("data" %in% slotNames(Z)){
      dt = Z@data
      return(noSpatCorr.test(Z=dt, X=coords, ...))
    }
  }
  stop("noSpatCorr.test: default method works only for Spatial***DataFrame objects.")
}


#' @describeIn noSpatCorr.test  Test for lack of spatial correlation
#' @method noSpatCorr.test matrix
#' @param R number of realizations of the Monte Carlo test
#' @param maxlag0 maximum lag distance to consider in the short range covariance
#' @param minlagInf minimum lag distance to consider in the long range covariance
#' @export
noSpatCorr.test.matrix <- function(Z, X, R = 299,
                       maxlag0=0.1*max(as.matrix(dist(X))), 
                       minlagInf=0.25*max(as.matrix(dist(X))),
                       ...){
  X = as.matrix(X)
  nG = ncol(X)
  n = nrow(X)
  X = unclass(idt(X))
    mystat = function(x, idx=1:n){
    X = x[,1:nG]
    Z = as.matrix(x[idx,-(1:nG)])
    ij = expand.grid(i=1:nrow(X), j=1:nrow(X))
    xd = (X[ij$i,]-X[ij$j,])^2 %*% rep(1, ncol(X))
    tk0 = xd<=(maxlag0^2)
    tk1 = xd>=(minlagInf^2)
    C0 = t(Z[ij[tk0,"i"],]-Z[ij[tk0,"j"],]) %*% (Z[ij[tk0,"i"],]-Z[ij[tk0,"j"],]) / (2*sum(tk0))
    C1 = t(Z[ij[tk1,"i"],]-Z[ij[tk1,"j"],]) %*% (Z[ij[tk1,"i"],]-Z[ij[tk1,"j"],]) / (2*sum(tk1))
    #ev = Re(eigen(solve(C0, C1), only.values = T)[[1]])
    #max(ev)/min(ev) + (max(ev)-1)^2
    C1sqrt = gsi.powM(C1, alpha=-0.5)
    ev = eigen(C1sqrt %*% C0 %*% C1sqrt, only.values = TRUE)[[1]]
    ev[length(ev)]
  }
  erg = boot::boot(cbind(X, Z), mystat, R=R, sim="permutation", ...)
  out = list(statistic = erg$t0, p.value = mean(erg$t<erg$t0), 
             method="permutation penalized eigenvalue ratio",
             maxlag0=maxlag0, minlagInf=minlagInf)
  class(out) = "htest"
  return(out)
}




# library("compositions")
# compileme=TRUE
# source("./R/geostats.R")
# turningbands = gsi.TurningBands
# CondSim = gsi.CondTurningBands
# # grid
# r=10
# x = seq(from=-4*r, to=+4*r, by=r/50)
# y = x
# xyGrid = expand.grid(X=x, Y=y)
# V = ilrBase(D=3, method = "balanced")
# #S1 = diag(c(0.5,1,2))
# S = diag(c(4,1)) #,5,1))
# S[1,2]<- S[2,1] <- 1
# M = r*diag(2)
# S = clrvar2ilr(ilrvar2clr(S, V=V))
# 
# #pdf("hello.pdf", width=7, height=7, fam="Times")
# #for(j in 1:3){
# vgram <- setCgram(type=vg.Exp,
#                   nugget=0*S,
#                   sill=S, #(nstru, nvar, nvar)
#                   anisRanges=M     # these are "classical" ranges (i.e. distances)
# )
# 
# plot(vgram) ### looks like we have interpreted M inversely in the plot!!!
# plot(vgram, vdir.up=c(0,1), col=2)
# plot(vgram, vdir.lo=c(1,0), vdir.up=c(0,1), col=2, lwd=3)
# par(plot(vgram, cov=FALSE))
# 
# 
# 
# rs = turningbands(as.matrix(xyGrid),vgram=vgram,nBands=400,nsim=1)
# 
# tk = sample(1:nrow(xyGrid), 1000)
# Xin = as.matrix(xyGrid[tk,])
# plot(Xin)
# Zin = as.matrix(rs[tk,,])
# dim(Zin)
# 
# vg = gsi.EVario2D(X = Xin, Z=Zin, azimuthNr = 1)
# 
# vg2 = gsi.EVario2D(X = Xin, Z=Zin, azimuthNr = 2)
# 
# plot(vg)
# plot(vg2, pch=19, col=2:3, lty=1)
# 
# 
# # Xin = matriz de coordenadas de input
# # Zin = matriz de valores de variables: cbind(z,z)
# CondSim(Xout=as.matrix(xyGrid),Xin,Zin,vgram=vgram,nbands=400,nsim=300) 
# noCorr.test(Xin, Zin, R=99)
# noCorr.test(Xin[sample.int(nrow(Xin), replace=TRUE),], Zin, R=99)  




#### diagonalisation measures
#' Compute diagonalisation measures
#' 
#' Compute one or more diagonalisation measures out of an empirical multivariate variogram.
#'
#' @param vgemp the empirical variogram to qualify
#' @param ... ignored
#'
#' @return an object of a similar nature to `vgemp`, but where the desired quantities are
#' reported for each lag. This can then be plotted or averages be computed.
#' @export
#' 
#' @details The three measures provided are 
#' \describe{
#' \item{absolute deviation from diagonality ("add")}{defined as the sum of all off-diagonal elements
#' of the variogram, possibly squared ($p=2$ if `quadratic=TRUE` the default; otherwise $p=1$)}}
#' \deqn{
#' \zeta(h)=\sum_{k=1}^n\sum_{j\neq k}^n \gamma_{k,j}^p(h)
#' }
#' \describe{
#' \item{relative deviation from diagonality ("rdd")}{comparing the absolute sum of off-diagonal elements  
#' with the sum of the diagonal elements of the variogram, each possibly squared ($p=2$ if `quadratic=TRUE`; 
#' otherwise $p=1$ the default)}}
#' \deqn{
#' \tau(h)=\frac{\sum_{k=1}^n\sum_{j \neq k}^n |\gamma_{k,j}(h)|^p}{\sum_{k=1}^n|\gamma_{k,k}(h)|^p}
#' }
#' \describe{
#' \item{spatial diagonalisation efficiency ("sde")}{is the only one requiring `vgemp0`, because it compares
#' an initial state with a diagonalised state of the variogram system}}
#' \deqn{
#' \kappa(h)=1-
#' \frac{\sum_{k=1}^n\sum_{j \neq k}^n |\gamma_{k,j}(h)|^p}{\sum_{k=1}^n\sum_{j \neq k}^n |\gamma_{(0)k,j}(h)|^p }
#' }
#' 
#' The value of $p$ is controlled by the first value of `method`. That is, the results with `method=c("rdd", "add")` 
#' are not the same as those obtained with `method=c("add", "rdd")`, as in the first case $p=1$ and in the second case $p=2$.
#' 
#'
#' @examples
#' data("jura", package="gstat")
#' X = jura.pred[, 1:2]
#' Z = jura.pred[,-(1:6)]
#' gm1 = make.gmCompositionalGaussianSpatialModel(data=Z, coords=X, V="alr")
#' vg1 = variogram(as.gstat(gm1)) 
#' (r1 = spatialDecorrelation(vg1, method=c("add", "rdd")))
#' plot(r1)
#' mean(r1)
#' require("compositions")
#' pc = princomp(acomp(Z))
#' v = pc$loadings
#' colnames(v)=paste("pc", 1:ncol(v), sep="")
#' gm2 = make.gmCompositionalGaussianSpatialModel(data=Z, coords=X, V=v, prefix="pc")
#' vg2 = variogram(as.gstat(gm2)) 
#' (r2 = spatialDecorrelation(vg2, method=c("add", "rdd")))
#' plot(r2)
#' mean(r2)
#' (r21 = spatialDecorrelation(vg2, vg1, method="sde") )
#' plot(r21)
#' mean(r21)
spatialDecorrelation = function(vgemp, ...) UseMethod("spatialDecorrelation",vgemp)
  

#' @describeIn spatialDecorrelation Compute diagonalisation measures
#' @method spatialDecorrelation gstatVariogram
#' @export
#' @param vgemp0 optionally, a reference variogram (see below; necessary for `method="sde"`)
#' @param method which quantities are desired? one or more of c("rdd", "add", "sde")
#' @param quadratic should the quantities be computed for a variogram or for its square? see below
spatialDecorrelation.gstatVariogram <- function(vgemp, vgemp0=NULL, method="add",quadratic=method[1]!="rdd",...){
 if(length(method)==1){  
  if(method=="add"){
    if(quadratic) vgemp$gamma = vgemp$gamma^2
    aux = split(vgemp, vgemp$id)
    erg = aux[[1]]
    erg$gamma = 0*erg$gamma
    tk = grep(".", names(aux), fixed=TRUE)
    for(o in tk) erg$gamma = erg$gamma + 2*abs(aux[[o]]$gamma)
   }else if(method=="rdd"){
     if(quadratic) vgemp$gamma = vgemp$gamma^2
     aux = split(vgemp, vgemp$id)
     erg = aux[[1]]
     erg$gamma = 0*erg$gamma
     ergdenom = erg$gamma
     tk = grep(".", names(aux), fixed=TRUE)
     for(o in tk) erg$gamma = erg$gamma + 2*abs(aux[[o]]$gamma)
     for(o in (1:length(aux))[-tk]) ergdenom = ergdenom + abs(aux[[o]]$gamma)
     erg$gamma = erg$gamma/ergdenom 
  }else if(method=="sde"){
    if(is.null(vgemp0)) stop("spatialDecorrelation: for method 'sde' a variogram for an initial configuration vgemp0 is needed!")
    aux = spatialDecorrelation(vgemp, method="add", quadratic = quadratic)
    aux2 = spatialDecorrelation(vgemp0, method="add", quadratic = quadratic)
    erg = aux
    erg$gamma = 1-aux$gamma/aux2$gamma
  }else stop("spatialDecorrelation: unkonw method! choose one of 'add', 'rdd' or 'sde'")
  erg$id =  factor(erg$id[, drop=T], labels=method)
  class(erg) = unique(c("spatialDecorrelationMeasure", class(erg) ))
  return(erg)
 }else{
   aux = lapply(method , function(m){
     spatialDecorrelation(vgemp=vgemp, vgemp0=vgemp0, method=m,quadratic=quadratic,...)
   })
   erg = aux[[1]]
   for(i in 2:length(aux)) erg = rbind(erg, aux[[i]])
   class(erg) = unique(c("spatialDecorrelationMeasure", class(erg) ))
   return(erg)
 }  
} 

#' @describeIn spatialDecorrelation Compute diagonalisation measures
#' @method spatialDecorrelation logratioVariogram
#' @export
spatialDecorrelation.logratioVariogram <- function(vgemp,vgemp0=NULL, method="add",...){
  stop("spatialDecorrelation meaningless for logratioVariogram")
} 

#' @describeIn spatialDecorrelation Compute diagonalisation measures
#' @method spatialDecorrelation gmEVario
#' @export
spatialDecorrelation.gmEVario <- function(vgemp,vgemp0=NULL, method="add",...){
  stop("not yet implemented")
} 



#' Average measures of spatial decorrelation
#' 
#' Compute average measres of spatial decorrelation out of the result of a call 
#' to [spatialDecorrelation()]
#'
#' @param x the outcome of a call to [spatialDecorrelation()]
#' @param ... further arguments to mean
#'
#' @return a mean for each measure available on `x`, i.e. this can be a vector.
#' The output is named.
#' @export
#' @seealso [spatialDecorrelation()] for an example
mean.spatialDecorrelationMeasure = function(x, ...){
  sapply(split(x$gamma, as.factor(x$id)), mean, ...)
} 


