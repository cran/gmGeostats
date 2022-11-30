#### gstat easy/easier interface for multivariate data
# setOldClass("gstat") ### breaks down

#' Fit an LMC to an empirical variogram
#' 
#' Fit a linear model of coregionalisation to an empirical variogram
#'
#' @param v empirical variogram
#' @param ... further parameters
#' @export
#' @return Method fit_lmc.gstatVariogram is a wrapper around [gstat::fit.lmc()], that calls this function
#' and gives the resulting model its appropriate class (c("variogramModelList", "list")). 
#' Method fit_lmc.default returns the fitted lmc (this function currently uses gstat as a 
#' calculation machine, but this behavior can change in the future)
#' @aliases fit_lmc fit_lmc.default fit_lmc.logratioVariogramAnisotropy
#'
#' @examples
#' data("jura", package = "gstat")
#' X = jura.pred[,1:2]
#' Zc = jura.pred[,7:13]
#' gg = make.gmCompositionalGaussianSpatialModel(Zc, X, V="alr", formula = ~1)
#' vg = variogram(gg)
#' md = gstat::vgm(model="Sph", psill=1, nugget=1, range=1.5)
#' gg = fit_lmc(v=vg, g=gg, model=md)
#' variogramModelPlot(vg, model=gg)
fit_lmc  <- function(v, ...) UseMethod("fit_lmc", v)



#' @describeIn fit_lmc wrapper around gstat::fit.lmc method
#' @param g spatial data object, containing the original data
#' @param model LMC or variogram model to fit
#' @param fit.ranges logical, should ranges be modified? (default=FALSE)
#' @param fit.lmc logical, should the nugget and partial sill matrices be modified (default=TRUE) 
#' @param correct.diagonal positive value slightly larger than 1, for multiplying the direct variogram 
#' models and reduce the risk of numerically negative eigenvalues
#' @export
#' @method fit_lmc gstatVariogram
fit_lmc.gstatVariogram <- function(v, g, model, fit.ranges = FALSE, fit.lmc = !fit.ranges, correct.diagonal = 1.0, ...){
  res = gstat::fit.lmc(as.gstatVariogram(v, ...), as.gstat(g, ...), as.variogramModel(model, ...),
                       fit.ranges = fit.ranges, fit.lmc = fit.lmc, correct.diagonal=correct.diagonal)
  class(res$model) = c("variogramModelList", "list")
  return(res)
}


#' @describeIn fit_lmc flexible wrapper method for any class for which methods 
#' for [as.gstatVariogram()], [as.gstat()] and [as.variogramModel()] exist.
#' In the future there may be direct specialised implementations not depending on
#' package gstat.
#' @export
#' @method fit_lmc default
fit_lmc.default <- function(v, g, model,...){
  origclass = class(g)
  res = fit_lmc(as.gstatVariogram(v), as.gstat(g), as.variogramModel(model), ...)$model
  # activate in the future
  res = as(res, origclass)
  return(res)
}
  

#' @describeIn fit_lmc method for logratioVariogram wrapping compositions::fit.lmc.
#' In the future there may be direct specialised implementations, 
#' including anisotropy (not yet possible).
#' @export
#' @method fit_lmc logratioVariogram
fit_lmc.logratioVariogram <- function(v, g, model,...){
  res = compositions::fit.lmc(as.logratioVariogram(v), as.CompLinModCoReg(model), ...)
  return(res)
}




#' Convert a regionalized data container to gstat
#' 
#' Convert a regionalized data container to a "gstat" model object
#' 
#' @param object regionalized data container
#' @param ... accessory parameters (currently not used)
#'
#' @return A regionalized data container of class "gstat", 
#' eventually with variogram model included. See [gstat::gstat()] for more info.
#' @aliases as.gstat.default 
#' @export
#'
#' @examples
#' data("jura", package = "gstat")
#' X = jura.pred[,1:2]
#' Zc = jura.pred[,7:13]
#' gg = make.gmCompositionalGaussianSpatialModel(Zc, X, V="alr", formula = ~1)
#' as.gstat(gg)
as.gstat <- function(object, ...) UseMethod("as.gstat", object)

#' @describeIn as.gstat default does nothing
#' @method as.gstat default
as.gstat.default <- function(object, ...){
  return(object)
} 
  
setGeneric("as.gstat", as.gstat)


# packs a regionalized composition and their geographic coordinates into a
#   gstat object after an appropriate logratio representation
#   coords: geographic coordinates (it works sure with data.frame)
#   compo: an acomp object (NOT TRANSFORMED)
#   V: can be either the matrix PSI (of the notes) or the strings "clr", "ilr" or "alr"
#   nscore: should data be marginally transformed to normal scores?
#   formulaterm: term for the formula argument of gstat (to control between UK and OK/SK)
#   ...: further arguments to gstat (e.g. for controlling neighbourhood or specyfing a mean for SK)
compo2gstatLR = function(coords, compo, V=ilrBase(compo), 
                         lrvgLMC=NULL, nscore=FALSE, 
                         formulaterm = "~1", prefix=NULL, ...){
  
  # prepare constants
  V0 = V
  D = ncol(compo)
  o = gsi.produceV(V=V,D=D,orignames=colnames(compo),giveInv=FALSE, prefix=prefix)
  prefix = o$prefix
  V = o$V
  # compute data (in lrs or in normal scores), set variable names
  Zlr = compositions::idt(compo, V=V)  
  if(nscore){
    source("nscore.R") # load the nscore.R functions
    prefix= paste("NS",prefix,sep="")
    Zlr = sapply(1:ncol(V), function(i){
      rs = nscore(Zlr[,i])
      aux = rs$nscore
      attr(aux,"trn.table") = rs$trn.table  # this ensures that the backtransformation is stored in the object
      return(data.frame(aux))
    })
    Zlr = as.data.frame(Zlr)
  }
  if(is.null(colnames(Zlr))) colnames(Zlr) = paste(prefix, 1:(D-1), sep="")
  # create gstat object
  spatdescr = paste("~",c(paste(colnames(coords),collapse=" + ")), sep="")
  gg = NULL
  for(i in 1:(D-1)){
    id = colnames(Zlr)[i]
    frm = paste(id, formulaterm, sep="")
    gg = gstat::gstat(gg, id=id, formula = stats::as.formula(frm), locations=stats::as.formula(spatdescr), data=data.frame(coords, Zlr), ...)
  }
  # if a logratio LMC was provided, convert it to gstat variogramModelList 
  #     and attach it
  if(!is.null(lrvgLMC) & !nscore){
    # space for a future conversion of variation-variogram models to gstat-LR-variograms
    gg$model = as.variogramModel(lrvgLMC, V=V0, prefix=prefix)
  }
  # return
  return(gg)
}

## version for rmult
rmult2gstat = function(coords, data, V="cdt", 
                         vgLMC=NULL, nscore=FALSE, 
                         formulaterm = "~1", prefix=NULL, ...){
  
  P = ncol(data)
  if(nscore){
    source("nscore.R") # load the nscore.R functions
    prefix= paste("NS",prefix,sep="")
    Z = sapply(1:P, function(i){
      rs = nscore(data[,i])
      aux = rs$nscore
      attr(aux,"trn.table") = rs$trn.table  # this ensures that the backtransformation is stored in the object
      return(data.frame(aux))
    })
    Z = as.data.frame(Z)
  }else{
    Z = data
  }
  if(is.null(colnames(Z))) colnames(Z) = paste(prefix, 1:P, sep="")
  # create gstat object
  spatdescr = paste("~",c(paste(colnames(coords),collapse=" + ")), sep="")
  gg = NULL
  for(i in 1:P){
    id = colnames(Z)[i]
    frm = paste(id, formulaterm, sep="")
    gg = gstat::gstat(gg, id=id, formula = stats::as.formula(frm), locations=stats::as.formula(spatdescr), data=data.frame(coords, Z), ...)
  }
  # if a logratio LMC was provided, convert it to gstat variogramModelList 
  #     and attach it
  if(!is.null(vgLMC) & !nscore){
    # space for a future conversion of variation-variogram models to gstat-LR-variograms
    gg$model = as.variogramModel(vgLMC, prefix=prefix)
  }
  # return
  return(gg)
}



getGstatData = function(gg # gstat object
){
  return(gg$data[[1]]$data@data)
}


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
#' occurred. This is useful for constructing complex diagrams, by giving 
#' argument `closeplot=FALSE` and then adding layers 
#' of information. If you want to "freeze" your plot, either give `closeplot=TRUE` or 
#' embed your call in another call to \code{\link{par}}, e.g. \code{par(variogramModelPlot(...))}.
#' @export
#' @importFrom gstat vgm
#' @seealso `gstat::plot.gstatVariogram()`
#' @method variogramModelPlot gstatVariogram
#' @family variogramModelPlot
#' @examples
#' data("jura", package="gstat")
#' X = jura.pred[,1:2]
#' Zc = jura.pred[,7:13]
#' gg = make.gmCompositionalGaussianSpatialModel(Zc, X, V="alr", formula = ~1)
#' vg = variogram(gg)
#' md = gstat::vgm(model="Sph", psill=1, nugget=1, range=1.5)
#' gg = fit_lmc(v=vg, g=gg, model=md)
#' variogramModelPlot(vg, model=gg)
variogramModelPlot.gstatVariogram = 
  function(vg,  # gstatVariogram object
           model = NULL,   # gstat  or variogramModelList object containing a variogram model fitted to vg
           col = rev(rainbow(1+length(unique(vg$dir.hor)))),
           commonAxis = FALSE,
           newfig = TRUE,
           closeplot = TRUE,
           ...
  ){
    # capture graphical parameters
    dotlist = list(...) 
    # check class of gg object and extract variable names
    if(is.null(model)){
      noms = levels(vg$id)
      vrnames = sort(unique(unlist(strsplit(noms, split=".", fixed=TRUE))))
      model = lapply(noms, function(i) gstat::vgm(model="Sph", psill=0, range=1, nugget=0))
      names(model)=noms
    }else{
      if(is(model, "variogramModelList")){
        vrnames = names(model)
        vrnames = vrnames[-grep(".", vrnames, fixed=TRUE)]
      }else if(is(model, "gstat")){
        vrnames = names(model$data)
        model = model$model
      }else if(is(model,"variogramModel") ){
        vrnames = levels(model$id)[1]
      }else{
        ggtr = tryCatch(as.variogramModel(model))
        if(inherits(ggtr,"try-error")){
          stop("argument 'gg' must either be a variogramModel, a gstat object with a variogramModel, a variogramModelList object, or an object convertible to one") 
        }else{
          return(variogramModelPlot(vg=vg,  gg = ggtr, col = col, commonAxis = commonAxis, newfig = newfig, ...))
        }
      } 
    }
    d = length(vrnames)
    # plot empirical vario!
    opar = par()
    opar = par_remove_readonly(opar)
    
    if(closeplot) on.exit(par(opar))
    
    if(newfig) par(mfrow=c(d,d), mar=c(1,1,1,1)+0.5, oma=c(0,3,3,0))  
    for(i in 1:d){
      for(j in 1:d){
        if(i==j){
          noms = vrnames[i]
        }else{
          noms = paste(vrnames[c(i,j)], vrnames[c(j,i)], sep=".")
        }
        # take the part of the empirical variogram corresponding to this (pair of) variable(s)
        tk = vg$id %in% noms
        empvar = vg[tk,]
        # find relevant stats
        azimuths = unique(empvar$dir.hor)
        rgH = range(empvar$dist, na.rm=TRUE)
        rgVG = range(empvar$gamma, na.rm=TRUE)
        if(commonAxis){
          tk.row = grep( vrnames[i],vg$id)
          rgVG = range(vg[tk.row,"gamma"], na.rm=TRUE)
        } 
        # take the corresponding model
        modvar = model[which(names(model) %in% noms )][[1]]
        # predict the several directions
        myfun = function(azimuth){
          dir = c(sin(azimuth*pi/180), cos((azimuth*pi/180)),0)
          gstat::variogramLine(modvar, max(rgH), dir=dir)
        }
        preds = lapply(azimuths, myfun)
        rgVG = range(0,rgVG, sapply(preds, function(X)X[,2]))
        ldotlist = dotlist
        if(!("xlim" %in% names(ldotlist))){
          ldotlist$xlim = range(0,empvar$dist, na.rm=TRUE) 
        }
        if(!("ylim" %in% names(ldotlist))){
          ldotlist$ylim = range(ifelse(i==j,0,NA),rgVG, na.rm=TRUE) 
        }
        if(!("pch" %in% names(ldotlist))){
          ldotlist$pch = 19 
        }
        if(!("lty" %in% names(ldotlist))){
          ldotlist$lty = 2
        }
        if(!("lwd" %in% names(ldotlist))){
          ldotlist$lwd = 2
        }
        if(!("xlab" %in% names(ldotlist))){
          ldotlist$xlab = ""
        }
        if(!("ylab" %in% names(ldotlist))){
          ldotlist$ylab = ""
        }
        ldotlist$col=col[as.integer(as.factor(empvar$dir.hor))]
        ldotlist$x = empvar$dist 
        ldotlist$y = empvar$gamma
        ldotlist$type="p"
        do.call(plot, args=ldotlist)
        sapply(1:length(azimuths), function(k){
          intk = empvar$dir.hor==azimuths[k]
          lines(gamma~dist, empvar[intk,], col=col[k], lty=ldotlist$lty)
          lines(preds[[k]], col=col[k], lwd=ldotlist$lwd )
        })
        if(i==1) mtext(side=3, text = vrnames[j], line=2 )
        if(j==1) mtext(side=2, text = vrnames[i], line=2 )
      }  
    }
    invisible(opar)
  }




#### functions to change between LMC and empirical variograms from/to gstat -------


## as gmEVario
# @describeIn as.gmEVario
# @export
as.gmEVario.gstatVariogram = function(vgemp, ...) stop("not yet available")

## as.logratioVariogram (empirical) -------
# transforms a gstat empirical variogram into a logratioVariogram object (evtl. with anisotropy)
# @describeIn as.logratioVariogram
as.logratioVariogram.gstatVariogram = function(vgemp,  # gstatVariogram object, emprical logratio variogram
                                               V=NULL, # matrix or name of the logratio transformation used
                                               tol=1e-12, # tolerance for generalized inverse (eventually for clr case)
                                               orignames=NULL, # names of the original component
                                               symmetrize=FALSE, # do you want a whole circle of directions?
                                               ...
){
  # prepare dimensions, names and constants
  DD = length(levels(vgemp$id))
  D = (1+sqrt(1+8*DD))/2
  if(is.null(orignames)) orignames = paste("v", 1:D, sep="")
  if(length(orignames)!=D) stop("names provided not consistent with number of logratio variables. Did you forget the rest?")
  o = gsi.produceV(V=V, D=D, orignames=orignames, giveInv=TRUE)
  prefix = o$prefix
  W = o$W
  # separate each direction, if anisotropic vario (ATTENTION: 3D not yet supported)
  vg4dir = split(vgemp, vgemp$dir.hor)
  # function to build one logratioVariogram
  buildOneLogratioVariogram = function(vg){
    aux = split(vg, vg$id)
    N = nrow(aux[[1]])
    lrnames = unique(names(aux))
    lrnames = sort(lrnames[-grep(".", lrnames, fixed=TRUE)])
    h = array(aux[[1]][,2], dim=c(N, D, D), dimnames=list(NULL, orignames, orignames))
    n = array(aux[[1]][,1], dim=c(N, D, D), dimnames=list(NULL, orignames, orignames))
    v = array(0, dim=c(D-1, N, D-1), dimnames=list(lrnames, NULL, lrnames))
    vvns = strsplit(names(aux), ".", fixed=TRUE)
    for(i in 1:length(vvns)){
      vns = vvns[[i]]
      if(length(vns)==1){ # diagonal
        v[vns, ,vns] = aux[[i]]$gamma
      }else{#off-diagonal
        v[vns[1], ,vns[2]] = aux[[i]]$gamma
        v[vns[2], ,vns[1]] = aux[[i]]$gamma     
      }
    }
    d = D-1
    dim(v) = c(N*d,d)
    v = v %*% W
    dim(v) = c(d, N*D)
    v = t(W) %*% v
    dim(v) = c(D,N,D)
    # v = aperm(v, c(2,1,3))
    v = gmApply(v,2,clrvar2variation)
    dim(v) = c(D,D, N)
    dimnames(v) = list(orignames, orignames, NULL)
    v = aperm(v, c(3,1,2))
    erg = structure(list(vg=v, h=h, n=n), class="logratioVariogram")
  }
  # create all variograms on all directions
  res = sapply(vg4dir, buildOneLogratioVariogram)
  # if a whole circle of directions is desired...
  if(symmetrize){
    cn = colnames(res)
    res = cbind(res, res)
    colnames(res) = c(cn, 180+as.double(cn))
  }
  # prepare and return the  "logratioVariogramAnisotropy" object
  hh = res["h",1][[1]][,1,1]
  mndfh = mean(diff(hh))
  dists = (hh[-1]+hh[-length(hh)])/2
  attr(res, "dists") = c(0, dists, max(dists)+mndfh)
  class(res)=c("logratioVariogramAnisotropy", "logratioVariogram")
  return(res)
}



## as.gstatVariogram (empirical) -------

#' Represent an empirical variogram in "gstatVariogram" format 
#' 
#' Represent an empirical variogram in "gstatVariogram" format, from package "gstat"; see [gstat::variogram()]
#' for details.
#' 
#' @param vgemp empirical variogram of any kind
#' @param ... further parameters (for generic functionality)
#'
#' @return The function returns an object of class "gstatVariogram" containing the empirical variogram provided.
#' See `gstat::variogram()` for details.
#' @export
#'
#' @examples
#' data("jura", package = "gstat")
#' X = jura.pred[,1:2]
#' Zc = compositions::acomp(jura.pred[,7:13])
#' lrvg = gmGeostats::logratioVariogram(data=Zc, loc=X)
#' as.gstatVariogram(lrvg, V="alr")
as.gstatVariogram <- function(vgemp, ...) UseMethod("as.gstatVariogram", vgemp)

#' @describeIn as.gstatVariogram Represent an empirical variogram in "gstatVariogram" format
#' @method as.gstatVariogram default
as.gstatVariogram.default <- function(vgemp, ...) vgemp


#' @describeIn as.gstatVariogram Represent an empirical variogram in "gstatVariogram" format
#' @method as.gstatVariogram gmEVario
as.gstatVariogram.gmEVario <- function(vgemp,...) stop("not yet available")

#' @describeIn as.gstatVariogram Represent an empirical variogram in "gstatVariogram" format
#' @method as.gstatVariogram logratioVariogram
#' @export
#' @param V eventually, indicator of which logratio should be used (one of: a matrix of logcontrasts, or of the strings "ilr", "alr" or "clr")
#' @param dir.hor eventually, which horizontal direction is captured by the  variogram provided (seldom to be touched!)
#' @param dir.ver eventually, which vertical direction is captured by the  variogram provided (seldom to be touched!)
#' @param prefix prefix name to use for the variables created (seldom needed)
as.gstatVariogram.logratioVariogram = 
  function(vgemp,  # gstatVariogram object, emprical logratio variogram
           V=NULL, # matrix or name of the logratio transformation used
           dir.hor=0,
           dir.ver=0,
           prefix=NULL,
           ...){
    class(vgemp) = NULL
    orignames = dimnames(vgemp$vg)[[2]]
    D = dim(vgemp$vg)[2]
    o = gsi.produceV(V=V, D=D,  orignames = orignames, giveInv = F, prefix=prefix )
    V = o$V
    prefix = o$prefix
    newnames = paste(prefix, 1:(D-1), sep="")
    Vu = V %*% diag(1/sqrt(colSums(V^2)))
    hh = gmApply(vgemp$h, 1, clrvar2ilr, V=Vu^2)
    nn = gmApply(vgemp$n, 1, clrvar2ilr, V=Vu^2)
    vv = -0.5*apply(vgemp$vg, 1, clrvar2ilr, V=V)
    ids = outer(newnames, newnames, paste, sep=".")
    diag(ids) = newnames
    
    ordre = NULL
    for(i in nrow(ids):1){
      ordre = c(ordre, ids[1:i,i])
    }
    
    rownames(nn) = ids
    rownames(hh) = ids
    rownames(vv) = ids
    
    # data.frame: np, dist, gamma, dir.hor, dir.ver=0, id= factor
    erg = data.frame(np=c(t(nn[ordre,])), dist=c(t(hh[ordre,])), 
                     gamma=c(t(vv[ordre,])), 
                     dir.hor=rep(dir.hor, each=ncol(nn)), 
                     dir.ver=rep(dir.ver, each=ncol(nn)),
                     id=factor(rep(1:length(ordre), each=ncol(nn)), labels=ordre)
    )
    class(erg) = c("gstatVariogram", "data.frame")
    return(erg)
  }


#' @describeIn as.gstatVariogram Represent an empirical variogram in "gstatVariogram" format
#' @method as.gstatVariogram logratioVariogramAnisotropy
#' @export
as.gstatVariogram.logratioVariogramAnisotropy = 
  function(vgemp,  # gstatVariogram object, emprical logratio variogram
           V=NULL, # matrix or name of the logratio transformation used
           ...){
    class(vgemp) = NULL
    alphas = gsi.azimuth2angle(colnames(vgemp))
    erg = as.gstatVariogram.logratioVariogram(vgemp[,1], V=V, dir.hor=alphas[1],...)
    for(i in 2:length(alphas)){
      erg = rbind(erg, as.gstatVariogram.logratioVariogram(vgemp[,i], V=V, dir.hor=alphas[i],...))
    }
    class(erg) = c("gstatVariogram", "data.frame")
    return(erg)
  }




## as.variogramModel (LMC) -------
#' Convert an LMC variogram model to gstat format
#' 
#' Convert a linear model of coregionalisation to the format of package gstat. See [gstat::vgm()] for details.
#'
#' @param m variogram model
#' @param ... further arguments for generic functionality
#'
#' @return The LMC model specified in the format of package gstat, i.e. as the result
#' of using [gstat::vgm()]
#' @export
#' @importFrom gstat gstat
#'
#' @examples
#' data("jura", package = "gstat")
#' X = jura.pred[,1:2]
#' Zc = compositions::acomp(jura.pred[,7:13])
#' lrmd = compositions::CompLinModCoReg(formula=~nugget()+sph(1.5), comp=Zc)
#' as.variogramModel(lrmd, V="alr")
as.variogramModel <- function(m, ...)  UseMethod("as.variogramModel", m)

#' @describeIn as.variogramModel Convert an LMC variogram model to gstat format
#' @method as.variogramModel default
#' @export
as.variogramModel.default <- function(m, ...) m

#' @describeIn as.variogramModel Convert an LMC variogram model to gstat format
#' @method as.variogramModel gmCgram
#' @export
as.variogramModel.gmCgram = function(m, ...) stop("not yet available")


#' @describeIn as.variogramModel Convert an LMC variogram model to gstat format
#' @method as.variogramModel LMCAnisCompo
#' @export
#' @param V eventually, specification of the logratio representation to use 
#' for compositional data (one of: a matrix of log-contrasts to use, or else one of 
#' the strings "alr", "clr" or "ilr")
#' @param prefix optional, name prefix for the generated variables if a transformation is used
#' @param ensurePSD logical, should positive-definiteness be enforced? defaults to TRUE, which may 
#' produce several scary looking but mostly danger-free warnings
as.variogramModel.LMCAnisCompo <- function(m, V=NULL, prefix=NULL, ensurePSD=TRUE, ...){
  D = ncol(m[,1]$sill)
  d=D-1
  o = gsi.produceV(V=V, D=D, orignames = dimnames(m["sill",1][[1]])[[2]], giveInv = FALSE, prefix=prefix)
  V = o$V
  prefix = o$prefix
  # which combinations of variables do we have to consider?
  noms = paste(prefix, 1:d, sep="")
  combs = cbind(rep(1:d, times=d:1), matrix(1:d, ncol=d, nrow=d)[lower.tri(diag(d), diag=TRUE)] )
  # equivalence table of correlogram model names
  eqtabmodels = factor(c(nugget="Nug", exp="Exp", sph="Sph", gau="Gau"), levels=levels(vgm()[,1]))
  models = eqtabmodels[sapply(1:ncol(m), function(j) m[,j]$model)]
  # express all sill matrices in the desired logratio 
  for(j in 1:ncol(m)){
    aux = -0.5 * t(V) %*% m[,j]$sill  %*% V
    colnames(aux) <- rownames(aux) <- noms
    if(ensurePSD){
      e = eigen(aux)
      tol = e$values[1]*1e-12
      e$values[e$values<tol]=tol
      aux = e$vectors %*% diag(e$values) %*% t(e$vectors)
      warning("as.variogramModel.LMCAnisCompo: negative eigenvalues corrected")
    }
    m[,j]$sill = aux
  }
  # recode A in azimuth, range and range ratio
  # TODO: correct and extend to 3D (check consistency with `anis2D.par2A()` and `anish2Dist()`)
  anis = matrix(0, nrow=ncol(m), ncol=3)
  colnames(anis) = c("range", "ang1", "anis1")
  for(j in 1:ncol(m)){
    anis[j, "anis1"] <- sqrt(sum((m[,j]$A[,2])^2))
    anis[j, "range"] <- m[,j]$range/anis[j, "anis1"]
    anis[j, "ang1"] <- atan2(-m[,j]$A[2,1], m[,j]$A[1,1]) * 180/pi
  }
  anis = data.frame(anis, ang2=0, anis2=1, ang3=0, kappa=0.5*(models!="Nug"))
  # if(all(anis$anis1==1)) anis=NULL
  # function to build one case
  myfun = function(ij){
    i = combs[ij,1]
    j = combs[ij,2]
    sills = sapply(1:ncol(m), function(k) m[,k]$sill[i,j] )
    objecte = data.frame(model=models, psill=sills)
    if(!is.null(anis))  objecte = cbind(objecte, anis)
    # class(objecte) = c("variogramModel","data.frame" )
    nugrow = which(objecte$model=="Nug")
    nugget = ifelse(length(nugrow)>0, objecte[nugrow,"psill"],0)
    first = (1:nrow(objecte))[-nugrow][1]
    md = vgm(model=objecte[first,"model"], psill=objecte[first,"psill"], range=objecte[first,"range"],
             nugget=nugget, anis=unlist(objecte[first,c("ang1","ang2","ang3", "anis1", "anis2")]),
             kappa = objecte[first, "kappa"])
    if(length((1:nrow(objecte))[-nugrow])>1)
      for(kk in (1:nrow(objecte))[-nugrow][-1])
        md = vgm(add.to=md, model=objecte[kk,"model"], psill=objecte[kk,"psill"], range=objecte[kk,"range"],
                 anis=unlist(objecte[kk,c("ang1","ang2","ang3", "anis1", "anis2")]),
                 kappa = objecte[kk, "kappa"])
    return(md)
  }
  res = lapply(1:nrow(combs), myfun)
  namelist = sapply(1:nrow(combs), function(ij) ifelse(combs[ij,1]==combs[ij,2], noms[combs[ij,1]], paste( noms[combs[ij,1]],  noms[combs[ij,2]], sep=".")   ) )
  names(res) = namelist
  #rownames(res) = NULL
  class(res) = c("variogramModelList","list")
  return(res)
}




#' @describeIn as.variogramModel Convert an LMC variogram model to gstat format
#' @method as.variogramModel CompLinModCoReg
#' @export
#' @param V eventually, specification of the logratio representation to use 
#' for compositional data (one of: a matrix of log-contrasts to use, or else one of 
#' the strings "alr", "clr" or "ilr")
#' @param prefix optional, name prefix for the generated variables if a transformation is used
#' @param ensurePSD logical, should positive-definiteness be enforced? defaults to TRUE, which may 
#' produce several scary looking but mostly danger-free warnings
#' @importFrom compositions vgram.nugget vgram.cardsin vgram.exp vgram.gauss vgram.lin vgram.pow vgram.sph
as.variogramModel.CompLinModCoReg <- function(m, V="alr", prefix=NULL, ensurePSD=TRUE, ...){
  strucs = gsi.extractCompLMCstructures(m)
  D = ncol(strucs$sills[[1]])
  o = gsi.produceV(V=V, D=D, giveInv = FALSE, prefix=prefix)
  as.variogramModel(as.LMCAnisCompo(m, varnames=rownames(o$V)), V=o$V, prefix=o$prefix, ensurePSD=ensurePSD, ...)
}


## as.LMCAnisCompo (LMC) -------
#' @describeIn as.LMCAnisCompo Recast compositional variogram model to format LMCAnisCompo
#' @method as.LMCAnisCompo gstat
#' @export
as.LMCAnisCompo.gstat <- function(m,...) as.LMCAnisCompo(m$model, ...)

# @describeIn as.LMCAnisCompo Recast compositional variogram model to format LMCAnisCompo
gstat2LMCAnisCompo <- as.LMCAnisCompo.gstat


#' @describeIn as.LMCAnisCompo Recast compositional variogram model to format LMCAnisCompo
#' @method as.LMCAnisCompo variogramModelList
#' @export
as.LMCAnisCompo.variogramModelList <- 
  function(m, V=NULL, orignames=NULL, ...){
    # prepare constants
    DD = length(m)
    D = (1+sqrt(1+8*DD))/2
    if(is.null(orignames)) orignames = paste("v", 1:D, sep="")
    if(length(orignames)!=D) stop("names provided not consistent with number of logratio variables. Did you forget the rest?")
    o = gsi.produceV(V=V, D=D, orignames = orignames, giveInv = TRUE)
    W = o$W
    colnames(W) = orignames
    # extract the dimensions and lr-names
    lrnames = unique(names(m))
    lrnames = sort(lrnames[-grep(".", lrnames, fixed=TRUE)]) # consider only names without "."
    # extract the number of structures
    K = nrow(m[[1]])
    if(length(lrnames)!=(D-1))stop("dimensions of data implied from m and V do not fit!")
    # equivalence table of correlogram model names
    eqtabmodels = c(Nug="nugget", Exp="exp", Sph="sph", Gau="gau")
    # function that extracts one particular sill matrix from the object
    darstellung = function(m, i){
      sill = matrix(0, ncol=D-1, nrow=D-1)
      rownames(sill) <- colnames(sill) <- lrnames
      for(jk in 1:DD){
        vvns = strsplit(names(m), ".", fixed=TRUE)[[jk]]
        if(length(vvns)==1) vvns = rep(vvns, 2)
        sill[vvns[1], vvns[2]] <- sill[vvns[2], vvns[1]] <- m[[jk]]$psill[i]
      }
      return(sill)
    }
    setvariostructure = function(i){
      model = eqtabmodels[ m[[1]]$model[i] ]
      sill = clrvar2variation(t(W) %*% darstellung(m, i)  %*% W)
      range = m[[1]]$range[i] * m[[1]]$anis1[i]
      # A = anis2D_par2A(ratio=m[[1]]$anis1[i], angle=m[[1]]$ang1[i], inv=FALSE)
      A = anis_GSLIBpar2A(ratios=c(m[[1]]$anis1[i],m[[1]]$anis2[i]), 
                          angles=c(m[[1]]$ang1[i],m[[1]]$ang2[i],m[[1]]$ang3[i]) )
      rs = list(model=model, range=range, A=A, sill=sill)
      class(rs) = "variostructure"
      return(rs)
    }
    res = sapply(1:K, setvariostructure)
    class(res) = "LMCAnisCompo"
    return(res)
  }


#' @describeIn as.gmCgram Convert theoretical structural functions to gmCgram format
#' @method as.gmCgram variogramModelList
as.gmCgram.variogramModelList = function(m, ...){ 
    as.gmCgram(as.LMCAnisCompo(m, ...), ...)
  }



#' @describeIn as.gmCgram Convert theoretical structural functions to gmCgram format
#' @method as.gmCgram variogramModel
as.gmCgram.variogramModel = function(m, ...){
  # extract nugget
  isNugget = m$model=="Nug"
  if(any(isNugget)){
    nuggetValue = m[isNugget, "psill"]
    m = m[!isNugget,, drop=FALSE]
  }
  # extract model names
  modelName = gsi.validModels[paste("vg", m$model,sep=".")]
  # if any model name is not identified
  if(any(is.na(modelName))){
    stop("as.gmCgram.variogramModel: found an unidentified variogram model; check content of internal variable gsi.valiModels to see which models are permissible")
  }
  # otherwise, extract parametres
  tt = function(x) t(t(x))
  out = setCgram(type = modelName[1], nugget = tt(nuggetValue), sill = tt(m[1, "psill"]), anisRanges = 
             as.AnisotropyScaling(unlist(m[1, -(1:4)])), extraPar = m[1, "kappa"])
  if(nrow(m)>1){
    for(im in 1:nrow(m)){
      out = out + setCgram(type = modelName[im], sill = tt(m[im, "psill"]), anisRanges = 
                       as.AnisotropyScaling(unlist(m[im, -(1:4)])), extraPar = m[im, "kappa"])
      
    }
  }
  return(out)
} 


