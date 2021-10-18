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
vg = vgm(model="Exp", psill=1, nugget=0, range=1) # variogram model
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
    data=cdt(Zdtc), coords = Xdt, # always set V="clr" in such cases!
    formula = ~1 # for ordinary (co)kriging
    )

## ---- fig.width=7, fig.height=5-----------------------------------------------
theta.vg = gmGeostats::variogram(theta.gg)
plot(theta.vg)

## ---- fig.width=7, fig.height=5-----------------------------------------------
theta.md = gstat::vgm(model="Exp", range=1, psill=0.5)
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

