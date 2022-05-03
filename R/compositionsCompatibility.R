#### empirical logratio variogram -----------
#' Empirical logratio variogram calculation
#' 
#' Calculation of an empirical logratio variogram (aka variation-variogram)
#' 
#' @param data a data container (of class "gmSpatialModel") or a composition (of class "acomp")
#' @param ... extra arguments for generic functionality
#' @return The output of this function depends on the number of `azimuth` values provided: 
#' if it is one single value (or if you explicitly call `logratioVariogram.default`) the result is 
#' a list of class "logratioVariogram" with the following elements
#'  * vg:	A nbins x D x D array containing the logratio variograms
#'  * h:	A nbins x D x D array containing the mean distance the value is computed on.
#'  * n:	A nbins x D x D array containing the number of nonmissing pairs used for the corresponding value.
#' If `azimuth` is a vector of directions, then the result is a matrix-like list of  "logratioVariogram" objects.
#' Each element of the mother list (or column of the matrix) is the variogram of one direction.
#' The output has a compound class c("logratioVariogramAnisotropy", "logratioVariogram").
#' @export
#'
#' @examples
#' data("jura", package="gstat")
#' X = jura.pred[,1:2]
#' Zc = compositions::acomp(jura.pred[,7:13])
#' vg = logratioVariogram(data=Zc, loc=X)
#' class(vg)
#' summary(vg)
#' vg = logratioVariogram(data=Zc, loc=X, azimuth=c(0,90))
#' class(vg)
#' summary(vg)
logratioVariogram <- function(data, ...) UseMethod("logratioVariogram", data)


#' @describeIn logratioVariogram Empirical logratio variogram calculation
#' @method logratioVariogram default
#' @export
#' @param loc if `data` is a composition (or if `comp` is provided), spatial coordinates of the sampling locations 
#' @param comp an alias for `data`, provided for back-compatibility with compositions::logratioVariogram 
logratioVariogram.default <- function(data, loc, ..., comp=data){
  res = try(compositions::logratioVariogram(comp=acomp(comp), loc=loc, ...), silent=TRUE)
  if(class(res)!="try-error") return(res)
  res = compositions::logratioVariogram(data=acomp(comp), loc=loc, ...)
}
  


# @describeIn logratioVariogram Empirical logratio variogram calculation for compositional data
# @param azimuth which direction, or directions, are desired (in case of directional variogram)
# @param azimuth.tol which tolerance sould be used for directional variograms?
# @export
# @method logratioVariogram acomp
logratioVariogram_acomp <- function(data=comp, loc, ..., azimuth=0, azimuth.tol=180/length(azimuth), comp=NULL){
  if(length(azimuth)>1) return(gsi.logratioVariogramAnisotropy(data, as.matrix(loc), ..., azimuths=azimuth, azimuth.tol=azimuth.tol))
  logratioVariogram.default(data, loc, ...)    
}


#' @describeIn variogram_gmSpatialModel logratio variogram method (see [logratioVariogram()] for details)
#' @export
#' @param data the data container (see [gmSpatialModel-class] for details)
#' @param azimuth which direction, or directions, are desired (in case of directional variogram)
#' @param azimuth.tol which tolerance sould be used for directional variograms?
logratioVariogram_gmSpatialModel <- function(data, ..., azimuth=0, azimuth.tol=180/length(azimuth)){
  coords = sp::coordinates(data)
  comp = tryCatch(gsi.orig(data@data))
  if(inherits(comp,"try-catch")) stop("logratioVariogram.gmSpatialModel: provided data is not compositional!")
  if(length(azimuth)>1) return(gsi.logratioVariogramAnisotropy(comp, coords, ..., azimuths=azimuth, azimuth.tol=azimuth.tol))
  logratioVariogram.default(comp, loc=coords, ...)    
}

#' @describeIn logratioVariogram S4 generic for gmSpatialModel objects 
setGeneric("logratioVariogram", logratioVariogram)




#' Logratio variogram of a compositional data
#' 
#' gmGeostats reimplementation of the compositions::logratioVariogram function
#'
#' @inheritParams logratioVariogram
#' @inheritParams logratioVariogram_gmSpatialModel
#' @inheritParams logratioVariogram.default 
#'
#' @return a "logratioVariogram" object, or a "logratioVariogramAnisotropy" object
#' if you provide more than one `azimuth`. See [logratioVariogram()] for details and 
#' examples.
#' @export
setMethod("logratioVariogram", c(data="acomp"), def=function(data=comp, loc, ..., azimuth=0, azimuth.tol=180/length(azimuth), comp=NULL){
  if(length(azimuth)>1) return(gsi.logratioVariogramAnisotropy(data, as.matrix(loc), ..., azimuths=azimuth, azimuth.tol=azimuth.tol))
  logratioVariogram.default(data, loc, ...)    
})


#################################################
#### Function to compute an empirical anisotropic -----------
## variogram: BUGGY! anisotropy not right in "compositions"
logratioVariogramAnisotropy1 = function(comp, loc, # 
                                        maxdist=max(dist(loc))/2, nbins=10,  # h maximal value and nr of classes
                                        dists = seq(0,maxdist,length.out=nbins+1), # h classes to consider
                                        azimuths=(0:11)*30,   # azimuths to consider
                                        azimuth.tol=360/length(azimuths)){   # angular tolerance
  onevario = function(a){
    return(logratioVariogram(comp, loc, dists=dists, azimuth=a, azimuth.tol=azimuth.tol))
  }
  rs = sapply(azimuths, onevario)
  class(rs) = c("logratioVariogramAnisotropy", "logratioVariogram")
  colnames(rs) = azimuths
  attr(rs, "lags") = dists
  return(rs)
}
#################################################
# Alternative function: do it myself within R
# will only work with small datasets
gsi.logratioVariogramAnisotropy = function(comp, loc, # 
                                       maxdist=max(dist(loc))/2, nbins=10,  # h maximal value and nr of classes
                                       dists = seq(0,maxdist,length.out=nbins+1), # h classes to consider
                                       azimuths=(0:11)*30,   # azimuths to consider
                                       azimuth.tol=360/length(azimuths)){   # angular tolerance
  # calculate the distances between all points, in angle and in radius
  dr = as.matrix(dist(loc))
  diag(dr) = NA # do not consider comparison of a point to itself
  # da = outer(1:nrow(loc), 1:nrow(loc), function(i,j) atan2(y=loc[j,2]-loc[i,2], x=loc[j,1]-loc[i,1]) )
  da = outer(1:nrow(loc), 1:nrow(loc), function(i,j) atan2(x=loc[j,2]-loc[i,2], y=loc[j,1]-loc[i,1]) ) # because angles are given in azimuths
  # compute the angular distance to each class
  angdist = function(a,b) gmApply(cbind(abs(a-b), abs(a-b+2*pi), abs(a-b-2*pi)), 1, min) # construct an angular distance function
  azimuths_rad = azimuths * pi/180
  comparisonA = outer(c(da), azimuths_rad, "-" ) %% (2*pi)  # angular distances are modulo 2pi
  azimuth.tol = azimuth.tol * pi/180
  # place each radial distances in its class
  distmids = dists[-1] - diff(dists)/2
  comparisonR = sapply(c(dr), function(dd){
    aux = abs(dd-distmids)
    which(aux==min(aux, na.rm=TRUE))[1]
  })
  # compute all increments
  Zclr = clr(comp)
  ij = expand.grid(i=1:nrow(Zclr), j=1:nrow(Zclr))
  Zincr = compositions::rmult(Zclr[ij$i,] - Zclr[ij$j,])
  # accessory function: compute noncentred variance of the increment on each class
  momenttwo = function(x){
    mn = unclass(mean(compositions::rmult(x)))
    mn = outer(mn, mn)
    return(var(compositions::rmult(x))+mn)
  } 
  # accessory function: form equivalent results to logratioVariogram from the preceding results
  mylogratioVariogram = function(i_Zincr, i_dr, i_comparisonR){
    # stats:
    nn = table(i_comparisonR)
    hh = sapply(split(i_dr, i_comparisonR), mean)
    vg = sapply(split(compositions::rmult(i_Zincr), i_comparisonR), momenttwo)
    # convert from clr to variation
    D = ncol(i_Zincr)
    vg = gmApply(vg,2, function(x){
      dim(x) = c(D,D)
      return(compositions::clrvar2variation(x))
    })
    vg = t(vg)
    dim(vg) = c(length(nn), D, D)
    dimnames(vg) = list(NULL,colnames(i_Zincr),colnames(i_Zincr))
    # prepare results
    hh = rep(hh, D*D)
    dim(hh) = dim(vg)
    nn = rep(nn, D*D)
    dim(nn) = dim(vg)
    dimnames(hh) <- dimnames(nn) <- dimnames(vg)
    erg = list(vg=0.5*vg, h=hh, n=nn) # semi-variogram
    class(erg) = "logratioVariogram"
    return(erg)
    
  }
  onevario = function(i){
    tk = abs(comparisonA[,i]) <= (azimuth.tol)/2
    aux = mylogratioVariogram(compositions::rmult(Zincr[tk,]), c(dr)[tk, drop=FALSE], as.factor(comparisonR)[tk, drop=FALSE])
    class(aux) = "logratioVariogram"
    return(aux)
  }
  rs = sapply(1:length(azimuths), onevario)
  class(rs) = c("logratioVariogramAnisotropy","logratioVariogram")
  colnames(rs) = paste(azimuths,"N", sep="")
  attr(rs, "lags") = gsi.lagdists(dists)
  attr(rs, "directions") = gsi.azimuth(azimuths)
  attr(rs, "type") = "semivariogram"
  attr(rs, "env") = environment()
  return(rs)
}
#################################################



#### splitting functions -------
### ATTENTION: TO MAKE coherent with gmEVario
#' Subsetting of logratioVariogram objects
#'
#' @param x an object of class "logratioVariogram" or c("logratioVariogramAnisotropy", "logratioVariogram")
#' @param i index or indexes of lags to be kept (if positive) or removed (if negative)
#' @param j index or indexes of directions to be kept, only for objects of class c("logratioVariogramAnisotropy", "logratioVariogram")
#' @param ... extra arguments, ignored
#'
#' @return the selected variograms or lags, potentially of class "logratioVariogram" if only one direction is chosen
#' @export
#' @aliases `[.logratioVariogram`
#'
#' @examples
#' data("jura", package="gstat")
#' X = jura.pred[,1:2]
#' Zc = compositions::acomp(jura.pred[,7:9])
#' vg = logratioVariogram(data=Zc, loc=X)
#' vg[1]
#' vg = logratioVariogram(data=Zc, loc=X, azimuth=c(0,90))
#' vg[,1]
#' vg[1,1]
#' vg[1,]
"[.logratioVariogramAnisotropy"<- function(x,i=NULL,j=NULL,...){
  if(is.null(i)){
    y = unclass(x)
    y = y[,j]
    attr(y, "lags") = attr(x, "lags")
    class(y) = class(x)
    if(length(j)==1) class(y) = "logratioVariogram"
    return(y)
  }
  if(is.null(j)){
    y = gmApply(x, 2, function(xx){
      class(xx) = "logratioVariogram"
      xx[i] 
    })
    class(y) = class(x)
    return(y)
  }
  x[,j][i,]
}

"[.logratioVariogram"<- function(x,i,...){
  y = lapply(x, function(xx) xx[i,,])
  class(y) = class(x)
  dists = attr(y,"lags")
  if(!is.null(dists)){
    attr(y,"lags") = c(dists[1], dists[-1][i])
  }
  return(y)
}



#################################################
# Alternative function to logratioVariogram
# compute one pwlr-variogram (no anisotropy embedded or allowed)
logratioVariogram1 = function(comp, loc, maxdist=max(dist(loc))/2, 
                              nbins=10, 
                              dists = seq(0,maxdist,length.out=nbins+1)){
  # calculate the distances between all points, in radius
  dr = as.matrix(dist(loc))
  diag(dr) = NA # do not consider comparison of a point to itself
  # place each radial distances in its class
  distmids = dists[-1] - diff(dists)/2
  comparisonR = sapply(c(dr), function(dd){
    aux = abs(dd-distmids)
    which(aux==min(aux, na.rm=TRUE))[1]
  })
  nn = table(comparisonR)
  hh = sapply(split(c(dr), comparisonR), mean)
  # compute all increments
  Zclr = clr(comp)
  ij = expand.grid(i=1:nrow(Zclr), j=1:nrow(Zclr))
  Zincr = compositions::rmult(Zclr[ij$i,] - Zclr[ij$j,])
  # compute noncentred variance of the increment on each class
  momenttwo = function(x){
    mn = unclass(mean(compositions::rmult(x)))
    mn = outer(mn, mn)
    return(var(compositions::rmult(x))+mn)
  }
  vg = sapply(split(Zincr, comparisonR), momenttwo)
  # convert from clr to variation
  D = ncol(comp)
  vg = gmApply(vg,2, function(x){
    dim(x) = c(D,D)
    return(clrvar2variation(x))
  })
  vg = t(vg)
  dim(vg) = c(length(nn), D, D)
  dimnames(vg) = list(NULL,colnames(comp),colnames(comp))
  # prepare results
  hh = rep(hh, D*D)
  dim(hh) = dim(vg)
  nn = rep(nn, D*D)
  dim(nn) = dim(vg)
  erg = list(vg=vg, h=hh, n=nn)
  class(erg) = "logratioVariogram"
  return(erg)
}
#################################################

#################################################
#### image method for logratioVariogramAnisotropy ----------
## objects. It generates a matrix of images of
## anisotropy for each possible logratio

#' Plot variogram maps for anisotropic logratio variograms
#' 
#' Image method to obtain variogram maps for anisotropic logratio variograms
#' 
#' @param x object of class c("logratioVariogramAnisotropy", "logratioVariogram")
#' @param jointColor logical, should all variogram maps share the same color scale?
#' @param breaks breaks to use in the color scale
#' @param probs alternatively to explicit `breaks`, these probabilities allow to ask for some equally probable breaks
#' @param col either a color palette, or else a vector of colors to use, of length `length(breaks)-1`
#' @param ... additional arguments for generic functionality (currently ignored)
#'
#' @return This function is called for its effect of producing a figure. 
#' Additionally, the graphical parameters active *prior to calling this function* are returned invisibly.
#' @export
#'
#' @examples
#' \donttest{
#' data("jura", package="gstat")
#' X = jura.pred[,1:2]
#' Zc = compositions::acomp(jura.pred[,7:9])
#' vg = logratioVariogram(data=Zc, loc=X, azimuth=c(0:18)*10, 
#'                            azimuth.tol=22.5)
#' image(vg)
#' image(vg, jointColor=TRUE)
#' }
image.logratioVariogramAnisotropy = function(x, jointColor=FALSE, breaks=NULL,
                                             probs = seq(0,1,0.1), col = spectralcolors,
                                             ...){
  lrvg = x
  o = par()
  # construct the polar grid
  r = attr(lrvg, "lags")
  rlim = range(r)
  phi = gsi.azimuth2angle(colnames(lrvg))
  dphi = mean(diff(phi))
  phi = phi - dphi/2
  phi = c(phi,phi[length(phi)]+dphi)*pi/180
  phi2 = pi/2-phi
  # extract the dimension of the composition
  class(lrvg) = NULL
  D = length(dimnames(lrvg[1,1][[1]])[[2]])
  # extract common z levels
  aux = c(unlist(sapply(1:ncol(lrvg), function(i) lrvg["vg",i][[1]])))
  if(is.null(breaks))
    breaks = quantile(aux[aux!=0], probs=probs, na.rm=TRUE)
  # set the matrix of figures
  par(mfrow=c(D,D), mar=c(1,1,1,1)/3)
  for(j in 1:D){
    for(k in 1:D){
      plot(c(-1,1)*rlim[2], c(-1,1)*rlim[2],type="n", asp=1, xlab="", ylab="", bty="n", xaxt="n", yaxt="n", ann=FALSE, xaxs="i",yaxs="i")
      if(j==k){
        text(0,0,dimnames(lrvg[1,1][[1]])[[2]][j], cex=2)      
      }else{
        # extract the directional variogram for the logratio k vs j
        z = sapply(1:ncol(lrvg), function(i){
          lrvg["vg",i][[1]][,j,k]
        })
        z[is.nan(z)]=NA
        dim(z) = c(length(r)-1, length(phi)-1)
        # plot it
        if(jointColor){
          if(is.function(col)) col=col(length(breaks)-1)
          image.polargrid(r,phi2,z,add=TRUE, breaks=breaks, col=col)
        }else{
          if(is.null(breaks)){
            if(is.function(col)) col=col(length(breaks)-1)
            image.polargrid(r,phi2,z,add=TRUE, col=col)
          }else{
            if(is.function(col)) col=col(length(probs)-1)
            image.polargrid(r,phi2,z,add=TRUE, col=col, probs=probs)
          }
        }  
      }
    }}   
  invisible(o)
}
#################################################



#################################################
#### plot method for logratioVariogramAnisotropy ------------
## objects. It generates a matrix of lines of
## variograms for each possible logratio and several
## directions. Potentially can include a model to
## copmpare with; as well as transform the variogram
## to other representations (clr, alr, ilr)
#' Plot variogram lines of empirical directional logratio variograms
#' 
#' Plots an "logratioVariogramAnisotropy" object in a series of panels, with each direction 
#' represented as a broken line.
#'
#' @param x logratio variogram with anisotropy, i.e. object of class c("logratioVariogramAnisotropy", "logratioVariogram")
#' @param azimuths which directions do you want to plot? default: all directions available
#' @param col colors to be used for plotting
#' @param type type of representation, see [graphics::plot()]
#' @param V optionally, a matrix of logcontrasts, or else one of the following strings: "alr", "ilr" or "clr"; 
#' to produce a plot of the empirical variogram in the corresponding representation; default to variation-variograms
#' @param lty style of the lines, potentially different for each directions
#' @param pch symbols for the points, potentially different for each directions
#' @param model eventually, variogram model to plot on top of the empirical variogram
#' @param figsp spacing between the several panels, if desired
#' @param closeplot boolean, should the plotting par() be returned to the starting values? (defaults to TRUE; 
#' see `plot.gmCgram()` for details) 
#' @param ... additional graphical arguments, to be passed to the underlying [graphics::matplot()] 
#' function
#'
#' @return Nothing. The function is called to create a plot.
#' @export
#' @importFrom stats predict
#'
#' @examples
#' data("jura", package="gstat")
#' X = jura.pred[,1:2]
#' Zc = compositions::acomp(jura.pred[,7:9])
#' vg = logratioVariogram(data=Zc, loc=X, azimuth=c(0:3)*45, 
#'                          azimuth.tol=22.5)
#' plot(vg) 
plot.logratioVariogramAnisotropy = function(x, azimuths=colnames(x), 
                                            col=rev(rainbow(length(azimuths))), type="o", 
                                            V = NULL, lty=1, pch=1:length(azimuths),
                                            model = NULL, figsp=0, closeplot=TRUE, ...){
  lrvg = x
  # construct the polar grid
  r = attr(lrvg, "lags")
  rlim = range(r)
  phi = gsi.azimuth2angle(colnames(lrvg))
  #phi = as.double(colnames(lrvg))
  azimuths = gsi.azimuth2angle(azimuths)
  # extract the dimension of the composition
  class(lrvg) = NULL
  D = length(dimnames(lrvg[1,1][[1]])[[2]])
  # check which representation is desired
  partnames = dimnames(lrvg["vg",1][[1]])[[2]]
  if(!is.null(V)){
    prefix = "ilr"
    if(is.character(V)){
      if(V=="ilr"){
        V = ilrBase(D=D)
        colnames(V) <- paste(prefix, 1:ncol(V), sep="")
      }else if(V=="alr"){
        V = rbind(diag(D-1), -1)
        prefix = "alr"
        colnames(V) <- paste(prefix, partnames[-D], sep="")
      }else if(V=="clr"){
        V = (diag(D)-matrix(1/D, ncol=D, nrow=D))#[, -D]
        prefix = "clr"
        colnames(V) <- paste(prefix, partnames, sep="")
      }
    }
    if(is.null(rownames(V))) rownames(V) <- partnames
  }else{
    V = diag(D)
    colnames(V) <- rownames(V) <- partnames
    prefix="variation"
  } 
  varnames = colnames(V)
  Dv = ncol(V)
  # if a model is provided, evaluate it
  if(!is.null(model)){
    rseq= seq(from=0, to=rlim[2]*1.03, length.out=200)
    myfun = function(azimuth){
      a = azimuth * pi/180
      vgd = predict(model, newdata=outer(rseq, c(sin(a), cos(a))) )
      if(prefix!="variation"){
        VvgdV = gmApply(vgd, 1, function(S) -0.5*t(V) %*% S %*% V)
        dim(VvgdV) <- c(Dv,Dv,length(rseq))
        vgd = aperm(VvgdV, c(3,1,2))
      } 
      dimnames(vgd)=list(NULL, colnames(V), colnames(V))
      return(vgd)
    }
    vgdAll = lapply(azimuths, myfun)
  }
  # transform the whole empirical variogram into the desired representation
  if(prefix!="variation"){
    for(i in 1:ncol(lrvg)){
      VvgdV = gmApply(lrvg[,i]$vg, 1, function(S) -0.5*t(V) %*% S %*% V)
      dim(VvgdV) <- c(Dv,Dv,dim(lrvg[,i]$vg)[1])
      vgd <- aperm(VvgdV, c(3,1,2))
      dimnames(vgd)=list(NULL, colnames(V), colnames(V))
      lrvg[,i]$vg <- vgd
    }
  } 
  
  # set the matrix of figures
  opar = par()
  opar = par_remove_readonly(opar)
  
  if(closeplot) on.exit(par(opar))
  
  
  par(mfrow=c(Dv,Dv), mar=figsp*c(1,1,1,1)/5, oma=c(3,3,3,3)+c(0,1,1,0)*ifelse(prefix=="variation",0,1), xaxs="i", yaxs="i")
  for(j in 1:Dv){
    vglim = range(sapply(1:ncol(lrvg), function(i){
      lrvg["vg",i][[1]][,j,]
    }), na.rm=TRUE)
    for(k in 1:Dv){
      if((j==k)&(prefix=="variation")){
        plot(c(0,1)*rlim[2], c(0,1)*vglim[2],type="n", asp=1, xlab="", ylab="", bty="n", xaxt="n", yaxt="n", ann=FALSE, xaxs="i",yaxs="i")
        text(0.5*rlim[2],0.5*vglim[2],dimnames(lrvg[1,1][[1]])[[2]][j], cex=2)      
      }else{
        # extract the directional variograms
        z = sapply(1:ncol(lrvg), function(i){
          lrvg["vg",i][[1]][,j,k]
        })
        z[is.nan(z)]=NA
        dim(z) = c(length(r)-1, length(phi))
        h = sapply(1:ncol(lrvg), function(i){
          lrvg["h",i][[1]][,j,k]
        })
        h[is.nan(h)]=NA
        dim(h) = c(length(r)-1, length(phi))
        # plot it
        tk = phi %in% azimuths
        if(!is.null(model)){
          vgloc = range(c(z[,tk],sapply(1:length(azimuths), function(iphi) vgdAll[[iphi]][,j,k])), na.rm=TRUE)
        }else{  
          vgloc = range(z[,tk], na.rm=TRUE)
        }
        ylm = range(0, ifelse(k>j,vglim[1],vgloc[1]), ifelse(k>j,vglim[2],vgloc[2]))
        matplot(h[,tk], z[,tk], col=col, type=type, xlim=c(0,1)*rlim[2], 
                ylim=ylm, xaxt="n", yaxt="n", lty=lty, pch=pch, ... )
        if(!is.null(model)){# add the model
          sapply(1:length(azimuths), function(iphi){
            lines(rseq, vgdAll[[iphi]][,j,k], lwd=2, col=col[iphi])
          })
        }
      }
      if(j %in% c(1,Dv) & j!=k) axis(side=ifelse(j==1,3,1))
      if(k %in% c(1,Dv) & j!=k) axis(side=ifelse(k==1,2,4))
      if(k==Dv & j!=k) axis(side=4) 
      #if(prefix=="variation" & (k-j)==1){
      #  axis(side=1)
      #  axis(side=2)
      #}
      #if(prefix=="variation" & (k-j)==(-1)){
      #  axis(side=3)
      #  axis(side=4)
      #}
      if(prefix!="variation" & j==1) mtext(side=3, text = varnames[k], line=2 )
      if(prefix!="variation" & k==1) mtext(side=2, text = varnames[j], line=2 )
    }}   
}

#' @describeIn fit_lmc method for logratioVariogram with anisotropry
#' @method fit_lmc logratioVariogramAnisotropy
#' @export
fit_lmc.logratioVariogramAnisotropy <- function(v, g, model, ...){
  warning("fit_lmc.logratioVariogramAnisotropy: attempting a workaround via gstat; expect changes in the future!")
  fit_lmc.gstatVariogram(v, g, model, ...)
}


#################################################


#################################################
## anisotropic variogram models and accesory functions 
# h : matrix of nx3 (or nx2) lag vectors
# A : isotropizing matrix
anish2Dist <- function (h, A = NULL){
  if (is.data.frame(h)) 
    h <- as.matrix(h)
  if(!is.null(A))
    if(is.matrix(A))
      if(nrow(A)>=ncol(h))
        h <- h %*% A[1:ncol(h),]
      if (is.matrix(h)) 
        h <- sqrt(rowSums(h^2))
      abs(h)
}

anvgram.exp = function (h, nugget = 0, sill = 1, range = 1, A = NULL, ...){
  "Exponential Variogram"
  s <- sill - nugget
  r <- range/log(10)
  h <- anish2Dist(h, A = A)
  ifelse(h > range * 1e-08, nugget + s * (1 - exp(-abs(h/r))), 0)
}



anvgram.gau = function (h, nugget = 0, sill = 1, range = 1, A = NULL, ...){
  "Gaussian Variogram"
  s <- sill - nugget
  r <- range/log(10)
  h <- anish2Dist(h, A = A)
  ifelse(h > range * 1e-08, nugget + s * (1 - exp(-(h/r)^2)), 0)
}

anvgram.sph = function (h, nugget = 0, sill = 1, range = 1, A = NULL, ...){
  "Spherical Variogram"
  h <- anish2Dist(h, A = A)
  s <- sill - nugget
  ifelse(h > range * 1e-08, nugget + s * ifelse(h > range, 
                                                1, 1.5 * h/range - 0.5 * (h/range)^3), 0)
}


anvgram.nugget = function (h, nugget = 0, sill = 0, range = 1, A = NULL, ...){
  "Nugget"
  anvgram.sph(h, nugget = sill, sill = 0, range = range, A = NULL, ...)
}




#' Create a anisotropic model for regionalized compositions  
#' 
#' Creates a (potentially anisotropic) variogram model for variation-variograms 
#'
#' @param Z compositional data set, used to derive the compositional dimension and colnames
#' @param models string (or vector of strings) specifying which reference model(s) to use 
#' @param azimuths typically a vector providing, for each model, the direction 
#' of maximal continuity (measured from North clockwise)
#' @param ranges typically a `G`-column matrix providing the minimal and maximal ranges, with one row per model
#' (with `G` specified below)
#' @param sillarray array of sills for each model. It can be null (to be estimated in the future). If specified, 
#' provide an appropriate `(x,x,K)`-array, where `K` is the number of models given, and `x` is explained below.
#' @param V depending on what do you give as sillarray, you may need to provide the matrix of logcontrasts, or
#' a string indicating whether the sills are represented in "alr" or "ilr"
#' @param tol tolerance for determination of positive definiteness
#' @param G one of c(1, 2, 3) identifying if we are in 1D, 2D or 3D cases
#'
#' @return an object of class "LMCAnisCompo" with all provided information appropriately structured, 
#' ready to be used for fitting, plotting or prediction. 
#' @export
#'
#' @examples
#' data("jura", package="gstat")
#' Zc = compositions::acomp(jura.pred[,7:9])
#' LMCAnisCompo(Zc, models=c("nugget", "sph"), azimuths=c(0,45))
LMCAnisCompo = function(Z, models=c("nugget", "sph", "sph"), 
                        azimuths=rep(0, length(models)), 
                        ranges=matrix(1, nrow=length(models), ncol=G), 
                        sillarray=NULL, V=NULL, tol=1e-12, G=2){
  # constants
  D = ncol(Z)
  d = D-1
  DD = D*(D+1)/2
  dd = d*D/2
  K = length(models)
  # check dimension of azimuths
  if(length(dim(azimuths))==0){
    dim(azimuths) = c(length(azimuths),1)
  }else if( !(ncol(azimuths) %in% c(1,3))){
    stop("LMCAnisCompo: `azimuths` should be a (column)vector in 2D, or a 3-column matrix in 3D")
  }
  
  # consider the sill matrices and try to interpret its shapes
  if(is.null(sillarray)){
    sillarray = rep(c(diag(D-1)), K)
    dim(sillarray) = c(D-1, D-1, K)
    sillarray = gmApply(sillarray,3,function(x) compositions::ilrvar2clr(x))
    dim(sillarray) = c(D, D, K)
  } 
  errormessage = "I cannot interpret sillarray; give either:\n
         a matrix of parameterPosdef(Clr)Mat's, or\n
         a (D,D,K) or (D-1, D-1,K) array, or\n
         a (D^2,K) or ((D-1)^2, K) matrix"
  if(!(length(sillarray) %in% c(D^2*K, d^2*K, dd*K, DD*K) ) )stop(errormessage)
  if(all(apply(sillarray, 3, compositions::is.clrvar))){
    darstellung = function(sa, i) compositions::clrvar2variation(sa[,,i])
  }else if(all(apply(sillarray, 3, compositions::is.ilrvar))){
    if(is.null(V)) stop("if you provide an ilr/alr representation, you must provide the matrix V as well")
    Vsvd = svd(V)
    W = with(Vsvd, v %*% diag(ifelse(d>tol, 1/d, 0)) %*% t(u) )
    darstellung = function(sa, i) compositions::clrvar2variation(t(W) %*% sa[,,i] %*% W)
  }else if(all(apply(sillarray, 3, compositions::is.variation))){
    darstellung = function(sa, i) sa[,,i]
  }else if(dim(sillarray)==c(DD,K)){
    darstellung = function(sa, i) clrvar2variation(parametricPosdefClrMat(sa[,i]))
  }else if(dim(sillarray)==c(dd,K)){
    if(is.null(V)) stop("if you provide an ilr/alr representation,you must provide the matrix V as well")
    Vsvd = svd(V)
    W = with(Vsvd, v %*% diag(ifelse(d>tol, 1/d, 0)) %*% t(u) )
    darstellung = function(sa, i) clrvar2variation(t(W) %*% parametricPosdefMat(sa[,i])%*% W)    
  }else{stop(errormessage)}    

  # consider isotropy
  if(length(dim(ranges))==0){
    ranges = matrix(rep(ranges, times=G), ncol=G)
  }
  
  # construct each structure
  setvariostructure = function(i){
    model = models[i]
    sill = darstellung(sillarray, i)
    range = ranges[i,2]
    # if(G==2)
    A = anis_GSLIBpar2A(ratios=ranges[i,-1]/ranges[i,1], angles=azimuths[i,] )
    # elseif(G==3) and so on...
    rs = list(model=model, range=range, A=A, sill=sill)
    class(rs) = "variostructure"
    return(rs)
  }
  res = sapply(1:K, setvariostructure)
  kk = grep(pattern = "nugget", x = models)
  if(length(kk)>0){
    res["A",kk]$A = 0*diag(ncol(res["A",kk]$A))
  }
  class(res) = "LMCAnisCompo"
  return(res)
}



# object must have the structure of a LMCAnisCompo
#    newdata is expected to contain lag vectors, i.e. N x 2(or 3) matrix or data.frame
#    function computes covariance values
#' Compute model variogram values
#'  
#' Evaluate the variogram model provided at some lag vectors
#' 
#' @param object variogram model
#' @param newdata matrix or data.frame of lag vectors
#' @param ... extra arguments for generic compatbility
#'
#' @return an array of dimension (nr of lags, D, D) with D the number of variables in the model.
#' @export
#' @include gmAnisotropy.R
#' @method predict LMCAnisCompo
#' @examples
#' data("jura", package="gstat")
#' Zc = compositions::acomp(jura.pred[,7:9])
#' lrmd = LMCAnisCompo(Zc, models=c("nugget", "sph"), azimuths=c(0,45))
#' predict(lrmd, outer(0:5, c(0,1)))
predict.LMCAnisCompo = function(object, newdata, ...){
  K = ncol(object)
  if(is.data.frame(newdata)) newdata = as.matrix(newdata)
  if(is.null(dim(newdata))) newdata = cbind(newdata, 0)
  if(ncol(newdata)==2) newdata = cbind(newdata,0)
  if(!(ncol(newdata) %in% c(2,3) )) stop("newdata should be 2D or 3D row-vectors")
  N = nrow(newdata)
  isnugget = which(unlist(object["model",])=="nugget")[1]
  ngt = object[,isnugget]$sill
  D = nrow(ngt)
  vvgg = rep(ngt,  N)
  dim(vvgg) = c(D, D, N)
  for(istruc in 1:K){
    fn = paste("anvgram", object[,istruc]$model, sep=".")
    g = eval(call(fn, h=newdata, range=object[,istruc]$range, sill=1, A=object[,istruc]$A))
    vvgg = vvgg + outer(object[,istruc]$sill, g)
  }
  vvgg = aperm(vvgg, c(3,1,2))
  return(vvgg)
}





#' Recast a model to the variogram model of package "compositions"
#' 
#' Recast a variogram model specified in any of the models of "gstat" or "gmGeostats" in
#' the format of [compositions::CompLinModCoReg()]
#'
#' @param v variogram model object to convert
#' @param ... further parameters for generic functionality
#'
#' @return The variogram model recast to "CompLinModCoReg" 
#' @export
#' @aliases as.CompLinModCoReg.CompLinModCoReg as.CompLinModCoReg.LMCAnisCompo
#' @importFrom compositions CompLinModCoReg 
as.CompLinModCoReg <- function(v, ...) UseMethod("as.CompLinModCoReg", v)

#' @method as.CompLinModCoReg CompLinModCoReg
#' @export
as.CompLinModCoReg.CompLinModCoReg <- function(v, ...)  v

#' @method as.CompLinModCoReg LMCAnisCompo
#' @export
as.CompLinModCoReg.LMCAnisCompo = function(v, ...){
  stop("as.CompLinModCoReg.LMCAnisCompo: not yet implemented")
}



# #################################################
# # function to fit one of the pwlr-variograms by clicking
# 
# clickfitOnePwlrVario = function(lrvg, i, j, variomodel, nclick=length(variomodel)+2,
#                                 exact =nclick==(length(variomodel)+2),
#                                 col=8){
#   if(class(lrvg)=="logratioVariogram"){
#     vg = lrvg$vg[,i,j]
#     h = lrvg$h[,i,j] 
#   }else if(class(lrvg)=="logratioVariogramAnisotropy"){
#     vg = sapply(1:ncol(lrvg), function(k) lrvg["vg",k][[1]][,i,j] )
#     h = sapply(1:ncol(lrvg), function(k) lrvg["h",k][[1]][,i,j] )
#   }
#   graphics::par(mfrow=c(1,1))
#   graphics::matplot(h, vg, type="l", col=col, lty=1, ylim=c(0, max(vg, na.rm=TRUE)), xlim=c(0, max(h, na.rm=TRUE)),
#           xaxs="i", yaxs="i")
#   xx = graphics::locator(nclick)
#   myvario = function(hh,par){
#     vvgg = rep(par[1], length(hh))
#     for(istruc in 1:length(variomodel)){
#       fn = paste("anvgram", variomodel[istruc], sep=".")
#       g = eval(call(fn, h=hh, range=par[istruc*2], sill=par[istruc*2+1]))
#       vvgg = vvgg + g
#     }
#     return(vvgg)
#   }
#   myfunerr = function(par){
#     sum((xx$y-myvario(xx$x, par))^2)
#   }
#   auxx = seq(from=min(xx$x[-1]), to=max(xx$x), length.out=length(variomodel))
#   auxy = seq(from=min(xx$y[-1]), to=max(xx$y), length.out=length(variomodel))-xx$y[1]
#   par0 = c(xx$y[1], c(cbind(auxx,auxy)))
#   max0 = rep(sapply(xx,max), times=length(variomodel))
#   max0 = c(max(xx$y), max0)
#   pp = stats::optim(par0, myfunerr, method="L-BFGS-B", lower=0*max0, upper=max0)
#   hdns = seq(from=0, to=max(h,na.rm=T), length.out=200)
#   lines(hdns, myvario(hdns, pp$par), col=2)
#   erg = list(models=variomodel, parameters=pp$par)
#   return(erg)
# }
# 
# 
# adjustOnePwlrVario = function(lrvg, i, j, variomodelpars,
#                               col=rainbow(ncol(lrvg)+3)){
#   requireNamespace("manipulate", quietly = TRUE)
#   if(is.null(col)) col=8
#   if(class(lrvg)=="logratioVariogram"){
#     vg = lrvg$vg[,i,j]
#     h = lrvg$h[,i,j] 
#     K = 1
#     # par(mfrow=c(1,1)) # not yet mature
#   }else if(class(lrvg)=="logratioVariogramAnisotropy"){
#     K = ncol(lrvg)
#     vg = sapply(1:K, function(k) lrvg["vg",k][[1]][,i,j] )
#     h = sapply(1:K, function(k) lrvg["h",k][[1]][,i,j] )
#     # par(mfrow=c(1,2)) # not yet mature
#   }
#   myvario = function(hh,par){
#     vvgg = rep(par[1], min(length(hh), nrow(hh)))
#     for(istruc in 1:length(variomodelpars$models)){
#       fn = paste("anvgram", variomodelpars$models[istruc], sep=".")
#       A = anis2D_par2A(par[istruc*4+2], par[istruc*4+3])
#       g = eval(call(fn, h=hh, range=par[istruc*4], sill=par[istruc*4+1], A=A))
#       # g = eval(call(fn, h=hh, range=par[istruc*4], sill=par[istruc*4+1]))
#       vvgg = vvgg + g
#     }
#     return(vvgg)
#   }
#   hdns = seq(from=0, to=max(h,na.rm=T), length.out=200)
#   par0 = variomodelpars$par
#   
#   doPlot = function(nugget, range1, sill1, ratio1, angle1, 
#                     range2, sill2, ratio2, angle2, 
#                     range3, sill3, ratio3, angle3){
#     requireNamespace("manipulate", quietly = TRUE)
#     par = c(nugget, range1, sill1, ratio1, angle1, range2, sill2, ratio2, angle2, range2, sill3, ratio3, angle3)
#     graphics::matplot(h, vg, type="l", col=col, lty=2, ylim=c(0, max(vg, na.rm=TRUE)), xlim=c(0, max(h, na.rm=TRUE)),
#             xaxs="i", yaxs="i" ) #main=sum((c(vg)-myvario(c(h), par))^2, na.rm=TRUE), 
#     for(k in 1:K){
#       if(K!=1){
#         aa = as.double(colnames(lrvg)[k]) * pi/180
#         vc = c(cos(aa), sin(aa))
#       }else{vc = c(1,0)}
#       graphics::lines(hdns, myvario(outer(hdns, vc), par), col=col[k], lwd=2, lty=1)
#     }
#     # imageOnePwlrAnis(lrvg, i, j) # not yet mature
#     manipulate::manipulatorSetState("OnePwlrPars", par)
#   }
#   mxvg = max(vg[is.finite(vg)], na.rm=T)
#   mxh =  max(h[is.finite(h)], na.rm=T)
#   manipulate::manipulate(
#     doPlot(nugget, range1, sill1, ratio1, angle1, 
#            range2, sill2, ratio2, angle2, 
#            range3, sill3, ratio3, angle3),
#     nugget = manipulate::slider(0, mxvg, initial=par0[1]),
#     range1 = manipulate::slider(0, mxh, initial=par0[2]),   
#     sill1 = manipulate::slider(0, mxvg, initial=par0[3]),
#     ratio1 = manipulate::slider(0, 1, initial=1),
#     angle1 = manipulate::slider(0, 180, initial=0),
#     range2 = manipulate::slider(0, mxh, initial=ifelse(length(par)>3, par0[4], 0)),   
#     sill2 = manipulate::slider(0, mxvg, initial=ifelse(length(par)>3, par0[5], 0)),   
#     ratio2 = manipulate::slider(0, 1, initial=1),
#     angle2 = manipulate::slider(0, 180, initial=0),
#     range3 = manipulate::slider(0, mxh, initial=ifelse(length(par)>5, par0[6], 0)),   
#     sill3 = manipulate::slider(0, mxvg, initial=ifelse(length(par)>5, par0[7], 0)),
#     ratio3 = manipulate::slider(0, 1, initial=1),
#     angle3 = manipulate::slider(0, 180, initial=0)
#   )  
# }
# 
# getOnePwlrPars = function(...){ manipulate::manipulatorGetState("OnePwlrPars") }
# 
# 
# adjustPwlrVario = function(lrvg, 
#                            cols=rainbow(ifelse(is.null(ncol(lrvg)),1, ncol(lrvg)))
# ){
#   par(mfrow=c(1,1))
#   requireNamespace("manipulate", quietly = TRUE)
#   if(class(lrvg)=="logratioVariogramAnisotropy"){
#     origdim = dim(lrvg["vg",1][[1]])
#     K = ncol(lrvg)
#     vg = sapply(1:K, function(k) lrvg["vg",k][[1]][,,] )
#     h = sapply(1:K, function(k) lrvg["h",k][[1]][,,] )
#     dim(h) <- dim(vg) <- c(origdim, K)
#   }
#   if(class(lrvg)=="logratioVariogram"){
#     vg = lrvg$vg[,,]
#     h = lrvg$h[,,] 
#     class(lrvg)="logratioVariogramAnisotropy"
#     lrvg = as.matrix(lrvg, ncol=1, nrow=3)
#     aux = h[,1,1] 
#     attr(lrvg, "lags") = c(0,aux)
#     colnames(lrvg) = "0"
#     dim(vg) <- dim(h) <- c(1,dim(h))
#   }
#   #myvario = function(hh,par){
#   #  vvgg = rep(par[1], length(hh))
#   #  for(istruc in 1:length(variomodelpars$models)){
#   #    fn = paste("anvgram", variomodelpars$models[istruc], sep=".")
#   #    g = eval(call(fn, h=hh, range=par[istruc*2], sill=par[istruc*2+1]))
#   #    vvgg = vvgg + g
#   #  }
#   #  return(vvgg)
#   #}
#   #hdns = seq(from=0, to=max(h,na.rm=T), length.out=200)
#   #par0 = variomodelpars$par
#   
#   myvario = function(hh,par){
#     vvgg = rep(par$pars[1], length(hh))
#     for(istruc in 1:length(par$models)){
#       fn = paste("anvgram", par$models[istruc], sep=".")
#       g = eval(call(fn, h=hh, range=par$pars[istruc*2], sill=par$pars[istruc*2+1]))
#       vvgg = vvgg + g
#     }
#     return(vvgg)
#   }
#   hdns = seq(from=0, to=max(h,na.rm=T), length.out=200)
#   
#   doPlot = function(i,j, nugget, model1, range1, sill1, model2, range2, sill2, model3, range3, sill3){
#     #manipulatorSetState("ij", c(i,j))
#     mm = c(model1, model2, model3)
#     nmod = mm !="none"
#     pp = c(nugget, range1, sill1, range2, sill2, range2, sill3)
#     par = list(models=mm[nmod] , pars=pp[c(TRUE,rep(nmod, each=2))])
#     #extras = manipulatorGetState
#     graphics::matplot(h[,i,j,], vg[,i,j,], type="l", col=cols, lty=1, ylim=c(0, max(vg, na.rm=TRUE)), xlim=c(0, max(h, na.rm=TRUE)),
#             xaxs="i", yaxs="i", main=sum((c(vg)-myvario(c(h), par))^2, na.rm=TRUE) )
#     graphics::lines(hdns, myvario(hdns, par), col=2)
#     D = dim(h)[2]
#     #if(loadIt){
#     #  TotalNugget = GetInitial("TotalNugget", D)
#     #  TotalSill1 =  GetInitial("TotalSill1", D)
#     #  TotalSill2 =  GetInitial("TotalSill2", D)
#     #  TotalSill3 =  GetInitial("TotalSill3", D)
#     #}  
#     #extras = list(TotalNugget, TotalSill1, TotalSill2, TotalSill3)
#     #if(saveIt){
#     #TotalNugget[j,i] <- TotalNugget[i,j] <- nugget
#     #TotalSill1[j,i] <- TotalSill1[i,j] <- sill1
#     #TotalSill2[j,i] <- TotalSill2[i,j] <- sill2
#     #TotalSill3[j,i] <- TotalSill3[i,j] <- sill3
#     #manipulatorSetState("TotalNugget") = parameterPosdefClrMat(variation2clrvar(TotalNugget))
#     #manipulatorSetState("TotalSill1") = parameterPosdefClrMat(variation2clrvar(TotalSill1))
#     #manipulatorSetState("TotalSill2") = parameterPosdefClrMat(variation2clrvar(TotalSill2))
#     #manipulatorSetState("TotalSill3") = parameterPosdefClrMat(variation2clrvar(TotalSill3))
#     #}
#   }
#   mxvg = max(vg[is.finite(vg)], na.rm=T)
#   mxh =  max(h[is.finite(h)], na.rm=T)
#   manipulate::manipulate(
#     doPlot(i,j, nugget, model1, range1, sill1, model2, range2, sill2, model3, range3, sill3),
#     i = manipulate::picker(as.list(1:(D-1))),
#     j = manipulate::picker(as.list((i+1):D)),
#     #loadIt = button("load"),
#     nugget = manipulate::slider(0, mxvg),
#     model1 = manipulate::picker("exp","sph"),
#     range1 = manipulate::slider(0, mxh),
#     sill1 = manipulate::slider(0, mxvg),
#     model2 = manipulate::picker("none", "exp","sph"),
#     range2 = manipulate::slider(0, mxh),   
#     sill2 = manipulate::slider(0, mxvg),
#     model3 = manipulate::picker("none", "exp","sph"),
#     range3 = manipulate::slider(0, mxh),   
#     sill3 = manipulate::slider(0, mxvg)  #,
#     # saveIt = button("save")
#   )  
# }
# 
# 
# SetInitial = function(statename){
#   ij = manipulate::manipulatorGetState("ij")
#   pp = manipulate::manipulatorGetState(statename)
#   if(!is.null(pp)){
#     aux = clrvar2variation(parametricPosdefClrMat(pp))
#     return(aux[ij[1], ij[2]])
#   }
# }
# 
# GetInitial = function(statename, D){
#   pp = manipulate::manipulatorGetState(statename)
#   if(!is.null(pp)){
#     aux = clrvar2variation(parametricPosdefClrMat(pp))
#     return(aux)
#   }else{
#     aux = matrix(0, nrow=D, ncol=D)
#     return(aux)    
#   }
# }







### manage the specification of basis confortably 
#' Create a matrix of logcontrasts and name prefix
#' 
#' Given a representation specification for compositions, this function
#' creates the matrix of logcontrasts and provides a suitable prefix name for naming variables. 
#'
#' @param V either a matrix of logcontrasts or, most commonly, one of "clr", "ilr" or "alr"
#' @param D the number of components of the composition represented 
#' @param orignames the names of the components
#' @param giveInv logical, is the inverse logcontrast matrix desired?
#' @param prefix the desired prefix name, if this is wished to be forced.
#'
#' @return A list with at least two elements
#'  1. `V` containing the final matrix of logcontrasts
#'  2. `prefix`  containing the final prefix for names of transformed variables
#'  3. `W` eventually, the (transposed, generalised) inverse of `V`, if `giveInv=TRUE`
#' @export
#'
#' @examples
#' gsi.produceV("alr", D=3)
#' gsi.produceV("ilr", D=3, orignames = c("Ca", "K", "Na"))
#' gsi.produceV("alr", D=3, orignames = c("Ca", "K", "Na"), giveInv = TRUE)
gsi.produceV = function(V=NULL, D=nrow(V), 
                        orignames=rownames(V), 
                        giveInv=FALSE, prefix=NULL){
  if(is.null(prefix)) prefix = "ilr"
  if(is.character(V)){
    if(V=="ilr"){
      V = ilrBase(D=D)
    }else if(V=="alr"){
      V = rbind(diag(D-1), -1)
      prefix = "alr"
    }else if(V=="clr"){
      V = (diag(D)-matrix(1/D, ncol=D, nrow=D))[, -D]
      prefix = "clr"
    }else if(V=="I"){
      V = diag(D)
      prefix=""
    }
  }
  if(!is.matrix(V)) stop("V must be an ilr-matrix or the strings 'ilr', 'alr', 'clr' or 'I'")

  if(is.null(rownames(V))){
    if(is.null(orignames)){
      rownames(V) = paste("z", 1:D, sep="")
    } else {
      rownames(V) = orignames
    }
  }
  
  if(is.null(colnames(V)))
    colnames(V) = paste(prefix, 1:ncol(V), sep="")
  
  if(giveInv){
    W = t(gsi.powM(V, alpha=-1))
    colnames(W) = rownames(V)
    rownames(W) = colnames(V)
    return(list(V=V, prefix=prefix, W=W))
  }
  return(list(V=V, prefix=prefix))
}





#' Quick plotting of empirical and theoretical logratio variograms
#' Quick and dirty plotting of empirical logratio variograms with or without their models 
#' @param vg empirical variogram or covariance function
#' @param model optional, theoretical variogram or covariance function
#' @param col colors to use for the several directional variograms
#' @param ... further parameters to `plot.logratioVariogramAnisotropy()`
#'
#' @return The function is primarily called for producing a plot. However, it 
#' invisibly returns the graphical parameters active before the call 
#' occurred. This is useful for constructing complex diagrams, by adding layers 
#' of info. If you want to "freeze" your plot, embed your call in another
#' call to \code{\link{par}}, e.g. \code{par(variogramModelPlot(...))}.
#' @export
#' @family variogramModelPlot
#' @method variogramModelPlot logratioVariogram
variogramModelPlot.logratioVariogram <- function(vg, model = NULL,   # gstat  or variogramModelList object containing a variogram model fitted to vg
                                        col = rev(rainbow(ndirections(vg))),
                                        ...){
  if(!is.null(model)){
    model = as.LMCAnisCompo(model)
  }
  vg = as.logratioVariogramAnisotropy(vg)
  plot(vg, col=col, model = model, ...)
} 



## as.gmEVario (empirical) -----
# @describeIn as.gmEVario
as.gmEVario.logratioVariogram = function(vgemp, ...) stop("not yet available")

# @describeIn as.gmEVario
as.gmEVario.logratioVariogramAnisotropy = function(vgemp, ...) stop("not yet available")




## as.logratioVariogram (empirical) -------

#' Recast empirical variogram to format logratioVariogram
#' 
#' Recast an empirical compositional variogram of any sort to a variation-variogram of class
#' "logratioVariogram".
#' 
#' @param vgemp empirical variogram
#' @param ... parameters for generic functionality
#'
#' @return the same model in the new format.
#' @export
#' @aliases as.logratioVariogram.logratioVariogram as.logratioVariogram.gmEVario
as.logratioVariogram <- function(vgemp,...) UseMethod("as.logratioVariogram",vgemp)

as.logratioVariogram.logratioVariogram <- function(vgemp, ...) vgemp

as.logratioVariogram.gmEVario  = function(vgemp, ...){ 
  stop("not yet available")
}




#' Convert empirical variogram to "logratioVariogramAnisotropy"
#' 
#' Convert an empirical variogram from any format to class "logratioVariogramAnisotropy" 
#'
#' @param vgemp an empirical variogram
#' @param ... further parameters
#'
#' @return The empirical variogram as a  "logratioVariogramAnisotropy" object 
#' @export 
as.logratioVariogramAnisotropy <- function(vgemp, ...) UseMethod("as.logratioVariogramAnisotropy", vgemp)

#' @describeIn as.logratioVariogramAnisotropy default method, making use of `as.logratioVariogram()`
#' @method as.logratioVariogramAnisotropy default
#' @export
as.logratioVariogramAnisotropy.default <- function(vgemp,...){
  rs = NextMethod(object=vgemp, ...)
  return(rs)
}

#' @describeIn as.logratioVariogramAnisotropy method for "logratioVariogram" class 
#' @method as.logratioVariogramAnisotropy logratioVariogram
#' @export
as.logratioVariogramAnisotropy.logratioVariogram <- function(vgemp,...){
  dim(vgemp) = c(3,1)
  attr(vgemp, "lags") = gsi.lagdists(apply(vgemp$h, 1, mean))
  attr(vgemp, "directions") = gsi.azimuth(0)
  colnames(vgemp) = "0N"
  class(vgemp) =  c("logratioVariogramAnisotropy", "logratioVariogram")
  attr(vgemp, "type") = "semivariogram"
  attr(vgemp, "env") = environment()
  return(vgemp)
}  
  

#' @describeIn as.logratioVariogramAnisotropy identity transformation
#' @method as.logratioVariogramAnisotropy logratioVariogramAnisotropy
#' @export
as.logratioVariogramAnisotropy.logratioVariogramAnisotropy <- function(vgemp,...){
  return(vgemp)
}




## as.LMCAnisCompo (LMC) -------

#' Recast compositional variogram model to format LMCAnisCompo
#' 
#' Recast a compositional variogram model of any sort to a variation-variogram model of class
#' "LMCAnisCompo".
#' 
#' @param m original variogram model
#' @param ... arguments for generic functionality
#' @return the variogram model recasted to class "LMCAnisCompo"
#' @aliases gstat2LMCAnisCompo
#' @export
as.LMCAnisCompo <- function(m, ...)  UseMethod("as.LMCAnisCompo",m)


#' @describeIn as.LMCAnisCompo Recast compositional variogram model to format LMCAnisCompo
#' @method as.LMCAnisCompo LMCAnisCompo
as.LMCAnisCompo.LMCAnisCompo <- function(m, ...) m


#' @describeIn as.LMCAnisCompo Recast compositional variogram model to format LMCAnisCompo
#' @method as.LMCAnisCompo gmCgram
#' @param V eventually, a specification of the way `m` is presently represented
#' @param orignames eventually, vector of names of the components, if `V` is provided and it does not have rownnames
as.LMCAnisCompo.gmCgram = function(m, V=NULL, orignames=rownames(V), ...){
  stop("not yet available")
} 


#' @describeIn as.LMCAnisCompo Recast a variogram model from package "compositions" to format LMCAnisCompo
#' @method as.LMCAnisCompo CompLinModCoReg
#' @param varnames a vector with the component names
#' @importFrom compositions parametricPosdefClrMat parametricRank1ClrMat
#' @export
as.LMCAnisCompo.CompLinModCoReg = function(m, varnames, ...){
  # extract parameters from the function
  strucs = gsi.extractCompLMCstructures(m)
  D = ncol(strucs$sills[[1]])
  if(missing(varnames)) varnames = paste("v", 1:D, sep="")
  if(D!=length(varnames)) stop("as.LMCAnisCompo: provides varnames does not have the appropriate length")
  # restructure sills as arrays
  sills = array(unlist(strucs$sills), dim = c(D,D, length(strucs$models)))
  # fake composition
  Z =  matrix(rep(1, D), nrow=1, dimnames = list("1", varnames))
  # output
  return( LMCAnisCompo(Z=Z, models=strucs$models, ranges=strucs$ranges, sillarray=sills) )
} 



#' Extract structures from a CompLinModCoReg model
#' 
#' Extract human-readable parameters from a CompLinModCoReg model
#'
#' @param m a CompLinModCoReg object
#'
#' @return A list with 3 elements ("models", "ranges" and "sills"), each in turn 
#' lists or vectors of the same length: "range" for nugget is NA; "models" are the string 
#' used in "compositions" for each model; "sills" is a list with the sills expressed
#' as variation matrices
gsi.extractCompLMCstructures = function(m){
  pars = compositions::vgmGetParameters(m)
  # split the parameter vector in the different nested structures 
  nm = names(pars)
  is = grep("r", nm)
  a = rep_len(0, length(pars))
  a[is]=1
  aa = 1+cumsum(a)
  strucs = split(pars, as.factor(aa))
  # extract the model names and the ranges
  aux = sapply(strucs, function(x) names(x)[1])
  models = gsub('[r,0-9]+', '', aux)
  models[models== "sPSD."] = "nugget"
  ranges = sapply(strucs, function(x) x[1])
  ranges[models=="nugget"] = NA
  # extract the sill matrices
  sills = split(pars[!a], as.factor( aa[!a]))
  types = sapply(sills, function(x) names(x)[1])
  types[grep("PSD", types)] = "parametricPosdefClrMat"
  types[grep("R", types)] = "parametricRank1ClrMat"
  M = lapply(1:length(sills), function(i)do.call(types[i], args = list(p=sills[[i]]))  )
  M = lapply(M, clrvar2variation)
  # done!
  return(list(models=models, ranges=ranges, sills=M))
}


## as.gmCgram (LMC) -------
# @describeIn as.gmCgram
as.gmCgram.LMCAnisCompo = function(m, V=attr(m,"contrasts"), 
                                   orignames=rownames(m["sill",1]$sill), ...){
  # produce constants
  utils::data("variogramModels")
  D = ncol(m["sill",1]$sill)
  d = D-1
  if(is.null(orignames)) orignames = paste("v", 1:D, sep="")
  if(length(orignames)!=D) stop("names provided not consistent with number of logratio variables. Did you forget the rest?")
  o = gsi.produceV(V=V, D=D, orignames = orignames, giveInv = FALSE)
  V = o$V
  prefix = o$prefix
  # prepare output
  output = list()
  # put structure models to output
  utils::data("variogramModels") # shortcut for all model constants
  eqtabmodels = c(gau=vg.Gauss, sph=vg.Sph, exp=vg.Exp)
  
  output$type = eqtabmodels[unlist(m["model",])]
  output$type = output$type[!is.na(output$type)]
  # put extra parameters (not active yet)
  output$data = 0*output$type
  # put nugget
  inugget = which(m["model",]=="nugget")
  if(length(inugget)!=0){
    output$nugget = -0.5*clrvar2ilr(m["sill",inugget]$sill, V=V)
  }else{
    output$nugget = 0*clrvar2ilr(m["sill",1]$sill, V=V)
  }
  # put anisotropy L matrix
  inonugget = which(m["model",]!="nugget")
  Ms = sapply(inonugget, function(i) as.AnisotropyRangeMatrix.AnisotropyScaling(m["A",i][[1]]) )
  dim(Ms) = c(dim(m["A",1][[1]]),length(inonugget))
  output$M = aperm(Ms, c(3,1,2))
  # put sills
  sills = sapply(inonugget, function(i) -0.5*clrvar2ilr(m["sill",i]$sill, V=V))
  dim(sills) = c(d,d,length(inonugget))
  output$sill = aperm(sills, c(3,1,2))
  class(output) = "gmCgram"
  # return output
  return(output)
}

# 
# # @describeIn as.variogramModel
#  EXISTS in gstatCompatibility.R
# as.variogramModel.CompLinModCoReg <- function(m, V=NULL, ...){
#   stop("not yet available")
#   # extract the substructures from m-variogram
#   rs = gsi.extractCompLMCstructures(m) # elements: "models", "ranges", "sills"
#   # construct the vgm template
#     # 1.- set the nugget (if needed) 
#   if(any(rs$models=="nugget")){
#     vg0 = vgm(psill=1, model="Nug")
#   }else{
#     vg0 = NULL
#   }
#    # 2.- add each structure
#   if(sum(rs$models!="nugget")>0){
#     for(i in which(rs$models!="nugget")){
#       vg0 = vgm(add.to = vg0, model=rs$models[i], range = rs$ranges[i], psill=1)
#     }
#   }
#   # cast the sills to the requested logratio coordinates
#   
#   # expand the vgm template to the new coordinates
#   
# } 


#' extract information about the original data, if available
#' 
#' originally implemented in package:compositions
#'
#' @param x an rmult object
#' @param y possibly another rmult object
#'
#' @return The original untransformed data (with a compositional class), resp
#' the TRANSPOSED INVERSE matrix, i.e. the one to recover the original components.
#' @aliases gsi.getV
gsi.orig <- function(x,y=NULL){
  a = attr(x,"orig")
  if(is.null(y)) return(a)
  b = attr(y,"orig")
  if(is.null(a)) return(b)
  return(a)
}
gsi.getV <- function(x,y=NULL){
  a = attr(x,"V")
  if(is.null(y)) return(a)
  b = attr(y,"V")
  if(is.null(a)) return(b)
  return(a)
}

