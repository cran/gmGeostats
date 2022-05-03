## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(compositions)
library(gstat)
library(gmGeostats)
library(magrittr)

## -----------------------------------------------------------------------------
circular = function(x, varname ="theta", conversion=pi/180){
  # output to be a (N, 1)-datamatrix 
  if(length(dim(x))!=2){ # case `x` is a vector
    y = t(t(x))
    colnames(y) = varname
  }else if(nrow(x)!=1){ # case `x` is a too large matrix
    y = x[, varname]
  }
  y = y * conversion
  class(y) = "circular"
  return(y)
}

## -----------------------------------------------------------------------------
cdt.circular = function(x, ...){ 
  z = cbind(sin(x), cos(x))
  colnames(z) = c("z1", "z2")
  return(rmult(z, orig=x))
}

## -----------------------------------------------------------------------------
cdtInv.circular = function(z, orig=compositions:::gsi.orig(x),...){ 
  #z = unclass(z)
  x = atan2(z[,2], z[,1])
  return(circular(x, varname=colnames(orig), conversion = 1))
}

## ----data_creation, fig.height=7, fig.width=7---------------------------------
## model setup
set.seed(333275)
xdt = data.frame(x=0, y=0, z=0) # one point is necesary
vg = vgm(model="Exp", psill=1, nugget=0, range=1, anis=c(30, 0.8)) # variogram model
gs = gstat(id="z", formula=z~1, locations = ~x+y , data=xdt, nmax=10, model=vg)
## sample point coordinates  
x <- runif(2000, min = 0, max = 10) # values between 0 and 10
Xdt = data.frame(x=x[1:1000], y=x[1001:2000])
# simulate random function
Z = predict(gs, newdata=Xdt, nsim=1)
# select columns
Zdt = Z[,3]
# define and plot data
Zdtc = circular(Zdt, varname = "theta", conversion = 1)
pairsmap(Zdtc, loc=Xdt)

## -----------------------------------------------------------------------------
theta.gg = 
  make.gmMultivariateGaussianSpatialModel(
    data=cdt(Zdtc), coords = Xdt, # always use cdt in such cases!
    formula = ~1 # for ordinary (co)kriging
    )

## ---- fig.width=7, fig.height=5-----------------------------------------------
theta.vg = gmGeostats::variogram(theta.gg)
plot(theta.vg)

## ---- fig.width=7, fig.height=5-----------------------------------------------
theta.md = gstat::vgm(model="Exp", range=1/3, psill=0.1) %>% 
  {gstat::vgm(add.to=., model="Sph", range=3, psill=0.1)}
theta.gs = fit_lmc(v=theta.vg, g = theta.gg, model = theta.md)
plot(theta.vg, model=theta.gs$model)

## -----------------------------------------------------------------------------
theta.gg = 
  make.gmMultivariateGaussianSpatialModel(
    data=cdt(Zdtc), coords = Xdt, # always set V="clr" in such cases!
    formula = ~1, # for ordinary (co)kriging
    model = theta.gs$model
    )

## -----------------------------------------------------------------------------
x <- seq(from=0, to=10, by=0.1)
xx = expand.grid(x,x)
colnames(xx) = colnames(Xdt)
ng = KrigingNeighbourhood(nmax = 20, omax=7, maxdist=1)
theta.prds = predict(theta.gg, newdata = xx, pars=ng)

## -----------------------------------------------------------------------------
theta.prds.grid = gsi.gstatCokriging2rmult(theta.prds)
theta.prds.back = backtransform(theta.prds.grid, as = cdt(Zdtc))
summary(theta.prds.back)

## ---- fig.width=6, fig.height=8.5---------------------------------------------
image_cokriged.circular = function(x, ...){
  class(x) = c("spatialGridRmult", "rmult")
  image_cokriged(x, breaks=40, col=rainbow(10), ...)
}
image_cokriged(theta.prds.back, ivar="theta")

## -----------------------------------------------------------------------------
classesToAM("gmSpatialMethodParameters", includeSubclasses = TRUE)

## ---- fig.width=7, fig.height=5-----------------------------------------------
# compute covariogram!
theta.cvg = gmGeostats::variogram(theta.gg, covariogram=TRUE)
class(theta.cvg)
# structure of the gstatVariogram object
head(theta.cvg)
# values controlling the split in direct and cross-variograms
theta.vg$id

# function doing the recalculations
recompute_complex_cov = function(cv){
  # split the gstatVariogram structure in the individual vgrams
  aux = split(cv[, -6], cv$id)
  # ad-hoc function taking two vgrams and operating their gamma column
  sumGamma = function(x,y, alpha=1){
    x$gamma = x$gamma + alpha*y$gamma
    return(x)
  }
  # compute C^{Im} and C^{Re}
  aaxx = list(Im=sumGamma(aux$z1.z2, aux$z1.z2, -1), 
              Re=sumGamma(aux$z1, aux$z2)
              )
  # undo the split
  f = rep(c("Im", "Re"), each=nrow(aux$z1))
  res = unsplit(aaxx, f)
  res$id = as.factor(f)
  # restore the class and return 
  class(res) = class(cv)
  return(res)
}
# do the calculations!
theta.vg %>% recompute_complex_cov %>% plot

## -----------------------------------------------------------------------------
# compute variogram for the whole circle, i.e. until 360 deg
theta.cvg = gmGeostats::variogram(
  theta.gg, covariogram=TRUE, alpha=(0:11)*30)
# how are the directions structured?
theta.cvg[theta.cvg$id=="z1", "dir.hor"]

## ---- fig.width=7, fig.height=5-----------------------------------------------
# function doing the recalculations
recompute_complex_cov_anis = function(cv){
  # split the gstatVariogram structure in the individual vgrams
  aux = split(cv[, -6], cv$id)
  # ad-hoc function taking two vgrams and operating their gamma column
  sumGamma = function(x,y, alpha=1){
    y$gamma = x$gamma + alpha*y$gamma
    return(y) # this time return the second argument!
  }
  N = nrow(aux[[1]])/2 
  # compute C^{Im} and C^{Re}
  aaxx = list(Im=sumGamma(aux$z1.z2[-(1:N),], aux$z1.z2[(1:N),], -1), 
              Re=sumGamma(aux$z1[1:N,], aux$z2[1:N,])
              )
  # undo the split
  f = rep(c("Im", "Re"), each=N)
  res = unsplit(aaxx, f)
  res$id = as.factor(f)
  # restore the class and return 
  class(res) = class(cv)
  return(res)
}

theta.cvg %>% recompute_complex_cov_anis() %>% plot

## -----------------------------------------------------------------------------
circularCovariogram = function(g, dirs=4){
  # case dirs is only the number of desired directions
  if(length(dirs)==1){
    delta = 180/dirs
    dirs = delta*(0:(dirs-1))
  }
  # extend the nr of directions to the whole circle
  dirs = c(dirs, 180 + dirs)
  # compute the covariogram
  cvg = gmGeostats::variogram(g, covariogram=TRUE, alpha=dirs)
  # recast and set the new class label
  rcvg = recompute_complex_cov_anis(cvg)
  class(rcvg) = c("circularCovariogram", class(cvg))
  return(rcvg)
}

## -----------------------------------------------------------------------------
# which arguments has as.gstatVariogram?
args(as.gstatVariogram)
# new method:
as.gstatVariogram.circularCovariogram = function(vgemp, ...){
  class(vgemp) = class(vgemp)[-1] # drop the extra class marker
  return(vgemp) # return result
}
# export the ad-hoc S3 class to S4
setOldClass("circularCovariogram")
# declare the new class an "EmpiricalStructuralFunctionSpecification"
setIs("circularCovariogram", "EmpiricalStructuralFunctionSpecification")

# check that all went right:
theta.cvg = circularCovariogram(theta.gg, dirs = 6)
is(theta.cvg, "circularCovariogram")
is(theta.cvg, "gstatVariogram")
is(theta.cvg, "data.frame")
is(theta.cvg, "EmpiricalStructuralFunctionSpecification")

