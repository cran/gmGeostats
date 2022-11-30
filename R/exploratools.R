#############################################################################
#####     FUNCTIONS FOR SPATIAL COMPOSITIONAL EXPLORATORY ANALYSIS     ######
#############################################################################



#' Check presence of missings
#' check presence of missings in a data.frame
#' @param x data, of class data.frame
#' @param mc not used
#' @param ... not used
#'
#' @return A single true/false response about the presence of any missing value 
#' on the whole data set
#' @export
#' @importFrom compositions has.missings
#' @method has.missings data.frame
#'
#' @examples
#' library(compositions)
#' data(Windarling)
#' has.missings(Windarling)
#' head(Windarling)
#' Windarling[1,1] = NA
#' head(Windarling)
#' has.missings(Windarling)
has.missings.data.frame = function(x, mc=NULL, ...){
  if(is.null(x))
    return(FALSE)
  (!is.null(x)) && any(is.na(x))
}




#' Compositional maps, pairwise logratios
#' Matrix of maps showing different combinations of components of a composition, in pairwise logratios
#' @param loc matrix or data.frame of coordinates of the sample locations
#' @param comp composition observed at every location, can be a matrix, a data.frame or
#'  of one of the classes \code{compositions::acomp} or \code{compositions::aplus}
#' @param colscale set of colors to be used as colorscale (defauls to 10 colors between blue and red)
#' @param cexrange symbol size min and max values (default to 0.1 to 2)
#' @param scale function scaling the set of z-values of each map, defaults to \code{\link{rank}}
#' @param commonscale logical, should all plots share a common z-scale? defaults to FALSE
#' @param foregroundcolor color to be used for the border of the symbol
#' @param closeplot logical, should the plot be left open (FALSE) for further changes, or be frozen (TRUE)? 
#' defaults to TRUE
#'
#' @return The function is primarily called for producing a matrix of (D,D) plots of the D-part 
#' compositional samples, where at each plot we represent a map whose symbols are colored and 
#' sized according to a z-scale controlled by a different logratio. For each plot, this is the 
#' logratio of the row variable by the column variable. However, in case that `closeplot=FALSE`, 
#' this function returns 
#' invisibly the graphical parameters that were active prior to calling this function. This allows 
#' the user to add further stuff to the plots (mostly, using \code{par(mfg=c(i,j))} to plot on the
#' diagram (i,j)), or manually freeze the plot (by wrapping the call to \code{pwlrmap} on a call to 
#' \code{\link{par}}).
#' @export
#' @importFrom graphics plot text mtext
#'
#' @examples
#' 
#' data("Windarling")
#' coords = as.matrix(Windarling[,c("Easting","Northing")])
#' compo = Windarling[,c("Fe","Al2O3","SiO2", "Mn", "P")]
#' compo$Rest = 100-rowSums(compo)
#' compo = compositions::acomp(compo) 
#' # in quantiles (default, ranking controls color and size)
#' pwlrmap(coords, compo) 
#' \donttest{
#' # in logratios (I=identity)
#' pwlrmap(coords, compo, scale=I)
#' # in ratios (i.e., apply exp)
#' pwlrmap(coords, compo, scale=exp)  
#' # use only color, no change in symbol size
#' pwlrmap(coords, compo, cexrange=c(1,1)) 
#' # change all
#' pwlrmap(coords, compo, commonscale=TRUE, cexrange=c(1.2,1.2), 
#'                     colscale=rev(rainbow(40, start=0, end=4/6))) 
#' }
pwlrmap = function(loc,   # XY coordinates (matrix or data frame)
                   comp,  # composition (matrix, acomp, aplus or data.frame)
                   colscale=rev(rainbow(10, start=0, end=4/6)), # color scale (defauls to 10 colors between blue and red)
                   cexrange=c(0.1, 2),  # symbol size min and max values (default to 0.1 to 2)
                   scale=rank, # scaling FUNCTION (defaults to ranking, i.e. quantiles)
                   commonscale=FALSE,  # should all plots be generated with a common scale?
                   foregroundcolor="black",
                   closeplot = TRUE
                   ){
  # set of maps where the symbols are chosen according to each possible pwlr, in 
  # a scale given by the user
  opar = par()
  opar = par_remove_readonly(opar)
  
  if(closeplot) on.exit(par(opar))
  # dimensions
  D = ncol(comp)
  N = nrow(loc)
  # compute pwlrs
  ij = expand.grid(i=1:D, j=1:D) # pairs of indices
  mat = matrix(0, nrow=D, ncol=D*D)   # matrix of +1 numerator, -1 denominator
  for(k in 1:nrow(ij)){
    mat[ij[k,1],k]=-1
    mat[ij[k,2],k]=1
  }
  Zpwlr = as.matrix(log(unclass(comp))) %*% mat
  # scale the variables
  if(commonscale){
    sizes  = scale(c(unlist(Zpwlr)))
    dim(sizes) = dim(Zpwlr)
    # calculate common cutting levels
    commonbks = seq(from=min(sizes[is.finite(sizes)]), to=max(sizes[is.finite(sizes)]), length.out=length(colscale))
    dfbks = diff(commonbks)[1]
    commonbks = c(commonbks[1]-dfbks, commonbks+0.5*dfbks)
  }else{
    sizes = gmApply(Zpwlr,2,scale)    
  }
  # do the plot
  par(mfrow=c(D,D), mar=c(1,1,1,1)/5, oma=c(3,5,5,3))
  for(i in 1:D){
    for(j in 1:D){
      if(i==j){
        # diagonal plots, only labels
        plot(loc, type="n", xaxt=ifelse(i==D,"s","n"), yaxt=ifelse(j==1,"s","n") )
        text(mean(range(loc[,1])), mean(range(loc[,2])), colnames(comp)[i], cex=2)
      }else{
      # off-diagonal plots, maps  
      k = (i-1)*D+j        
      sz = sizes[,k]  
      # choose which breaks to use, the common ones or the particular ones
      if(commonscale){
        bks = commonbks
      }else{
       bks = seq(from=min(sz[is.finite(sz)]), to=max(sz[is.finite(sz)]), length.out=length(colscale))
         dfbks = diff(bks)[1]
       bks = c(bks[1]-dfbks, bks+0.5*dfbks)
      }
      # compute colors and sizes
      cols = colscale[cut(sz, bks)]
      sz = cexrange[1]+(cexrange[2]-cexrange[1])*(sz-bks[1])/c(range(bks) %*% c(-1,1))
      # actual map
      plot(loc, cex=sz, bg=cols, pch=21, asp=1, col=foregroundcolor,
            xaxt=ifelse(i==D,"s","n"),yaxt=ifelse(j==1,"s","n")
           )
      }
      # add axes on extra sides
      if(i==1) axis(side=3)
      if(j==D) axis(side=4)
    }
  }
  # add labels
  mtext("numerator", side=2, outer=TRUE, line=2.5, cex=1.25)
  mtext("denominator", side=3, outer=TRUE, line=2.5, cex=1.25)
  # return the old graphical parameters to freeze the plot
  invisible(opar)
}






#' Multiple maps
#' Matrix of maps showing different combinations of components of a composition, user defined
#' 
#'
#' @param data data to represent; admits many data containing objects 
#' (matrix, data.frame, classes from package \code{compositions}) as well 
#' as their \code{Spatial} counterparts (in which case, \code{loc} is not necessary)
#' @param ... extra arguments for generic functionality
#'
#' @return The function is primarily called for producing a matrix of plots of each component of a 
#' multivariate data set, such as a compositional data set. Each plot represents a map whose symbols 
#' are colored and sized according to a z-scale controlled by one of the variables of the data set. 
#' It can be used virtually with any geometry, any kind of data (compositional, positive, raw, etc). 
#' This function returns invisibly the graphical parameters that were active prior to calling this 
#' function. This allows the user to add further stuff to the plots (mostly, using \code{par(mfg=c(i,j))} 
#' to plot on the diagram (i,j)), or else freeze the plot (by wrapping the call to \code{pwlrmap} 
#' on a call to \code{\link{par}}).
#' @export
#' @importFrom grDevices rainbow
#' @importFrom graphics par plot
#'
#' @examples
#' data("Windarling")
#' library("compositions")
#' coords = as.matrix(Windarling[,c("Easting","Northing")])
#' compo = Windarling[,c("Fe","Al2O3","SiO2", "Mn", "P")]
#' compo$Rest = 100-rowSums(compo)
#' compo = acomp(compo)
#' pairsmap(data=clr(compo), loc=coords) # clr
#' pairsmap(data=compo, loc=coords) # closed data
pairsmap <- function(data, ...) UseMethod("pairsmap", data)


#' @describeIn pairsmap Multiple maps
#' @method pairsmap SpatialPointsDataFrame
pairsmap.SpatialPointsDataFrame <- function(data, ...){
  pairsmap.default(data@data, loc=sp::coordinates(data), ...)
}


#' @describeIn pairsmap Multiple maps
#' @method pairsmap default
#' @export
#' @param loc matrix or data.frame of coordinates of the sample locations
#' @param colscale set of colors to be used as colorscale (defauls to 10 colors between blue and red)
#' @param cexrange symbol size min and max values (default to 0.1 to 2)
#' @param scale function scaling the set of z-values of each map, defaults to \code{\link{rank}}
#' @param commonscale logical, should all plots share a common z-scale? defaults to FALSE
#' @param mfrow vector of two integers, giving the number of plots per row and per column of the 
#' matrix of plots to be created; defaults to something square-ish, somewhat wider than longer, and able to
#' contain all plots
#' @param foregroundcolor color to be used for the border of the symbol
#' @param closeplot logical, should the plot be left open (FALSE) for further changes, or be frozen (TRUE)? 
#' defaults to TRUE
pairsmap.default <- function(data,   # data to represent
                    loc,   # XY coordinates (matrix or data frame)
                   colscale=rev(rainbow(10, start=0, end=4/6)), # color scale (defauls to 10 colors between blue and red)
                   cexrange=c(0.1, 2),  # symbol size min and max values (default to 0.1 to 2)
                   scale=rank, # scaling FUNCTION (defaults to ranking, i.e. quantiles)
                   commonscale=FALSE,  # should all plots be generated with a common scale?
                   mfrow = c(floor(sqrt(ncol(data))), ceiling(ncol(data)/floor(sqrt(ncol(data))))),# automatic distribution of figs
                   foregroundcolor = "black", 
                   closeplot=TRUE,
                   ...
){
  opar = par()
  opar = par_remove_readonly(opar)
  if(closeplot) on.exit(par(opar))
  # dimensions
  D = ncol(data)
  N = nrow(loc)
  # scale the variables
  if(commonscale){
    sizes  = scale(c(unlist(data)))
    dim(sizes) = dim(data)
    # calculate common cutting levels
    commonbks = seq(from=min(sizes[is.finite(sizes)]), to=max(sizes[is.finite(sizes)]), length.out=length(colscale))
    dfbks = diff(commonbks)[1]
    commonbks = c(commonbks[1]-dfbks, commonbks+0.5*dfbks)
  }else{
    sizes = gmApply(data,2,scale)    
  }
  # do the plot
  par(mfrow=mfrow, mar=c(1,1,10,1)/5, oma=c(3,3,3,3))
  drow = mfrow[1]
  dcol = mfrow[2]
  for(i in 1:drow){
    for(j in 1:dcol){
        # off-diagonal plots, maps  
        k = (i-1)*dcol+j        
        if(k<=D){
         sz = sizes[,k]  
         # choose which breaks to use, the common ones or the particular ones
         if(commonscale){
           bks = commonbks
         }else{
            bks = seq(from=min(sz[is.finite(sz)]), to=max(sz[is.finite(sz)]), length.out=length(colscale))
           dfbks = diff(bks)[1]
           bks = c(bks[1]-dfbks, bks+0.5*dfbks)
         }
         # compute colors and sizes
         cols = colscale[cut(sz, bks)]
         sz = cexrange[1]+(cexrange[2]-cexrange[1])*(sz-bks[1])/c(range(bks) %*% c(-1,1))
         # actual map
         plot(loc, cex=sz, bg=cols, pch=21, asp=1, main=colnames(data)[k], col=foregroundcolor,
             xaxt=ifelse(i==drow,"s","n"),yaxt=ifelse(j==1,"s","n")
         )
        } 
       # add axes on extra sides
       # if(i==1) axis(side=3)
       if(j==dcol) axis(side=4)
    }
  }
  # return the old graphical parameters
  #par(mfrow=opar$mfrow, mar=opar$mar, oma=opar$oma)
  invisible(opar)
}



#' Spectral colors palette
#' based on the RColorBrewer::brewer.pal(11,"Spectral")
#'
#' @param n number of colors 
#'
#' @return a palette, i.e. a list of colors, from dark blue to dark red over pale yellow.
#' @export
#' @importFrom grDevices colorRampPalette
#' @import RColorBrewer
#'
#' @examples
#' (cls=spectralcolors(20))
spectralcolors <- function(n){
  cls = RColorBrewer::brewer.pal(min(11,n), "Spectral")
  if(n>11){
    cls = grDevices::colorRampPalette(cls)(n)
  }
  return(rev(cls))
}

# #### MOVED TO COMPOSITIONS #####
# Compositional panel 1D-density plot
# Panel minifunction for plotting histograms and kernel densities of the data 
#
# @param x ignored (here for compatibility with \code{\link{qqnorm.acomp}})
# @param y numeric vector of response values
# @param col color of the plot
# @param ... further parameters to plotting functions, currently ignored
# @param alpha alpha level at which the test should be done
#
# @return If given to the argument \code{panel} of a function such as \code{\link{qqnorm.acomp}}),
# this produces a matrix of plots where each panel contains a histogram and a kernel density 
# overdrawn. If the distribution of this data is accepted to be normal at the \code{alpha}-level
# by a \code{\link{shapiro.test}}, then the histogram is painted with the \code{col}or provided; 
# otherwise the histogram bars are empty, but the kernel density curve is colored.
#
# @examples
# data("Windarling")
# compo = as.matrix(Windarling[,c("Fe","Al2O3","SiO2", "Mn", "P")])
# qqnorm.acomp(compo, panel=vp.lrdensityplot, alpha=0.05) 
# #### MOVED TO COMPOSITIONS #####





# Panel function for 2D-density plots
# Panel minifunction for plotting 2D kernel densities of the data 
#
# @param x numeric vector of response values (axis X)
# @param y numeric vector of response values (axis Y)
# @param xaxt style of the X axis labelling (defaults to "n", none)
# @param yaxt style of the Y axis labelling (defaults to "n", none)
# @param grid logical, should a unit grid be drawn? defaults to TRUE
# @param legpos string, position of the correlation coefficient, defaults to "bottomright"
# @param add logical, should the output be added to an existing diagram? defaults to TRUE, as required for acting as a panel diagram
# @param colpalette color palette for the image (defaults to spectral colors: blue-yellow-red)
# @param ... further parameters to image
#
# @return If given to the argument \code{panel} of a function such as \code{\link{pairs}}),
# this produces a matrix of images where each panel contains a 2D kernel density map,
# using blue for low density regions and dark red for high density colors. 
# Regression lines (y~x) and correlation coefficients are added as well.
# @export
#
# @examples
# data("Windarling")
# compo = Windarling[,c("Fe","Al2O3","SiO2", "Mn", "P")]
# pairs(clr(compo), panel=vp.kde2dplot)
# vp.kde2dplot = 
#   function(x, y, xaxt="n", yaxt="n",
#            grid=TRUE,legpos="bottomright", add=TRUE, colpalette=spectralcolors,...){
#     aux = MASS::kde2d(x, y, n=50)
#     aux$z = sqrt(aux$z)
#     bks = hist(aux$z, plot=F, breaks=20)$breaks
#     cols = c(NA,colpalette(length(bks)-2))
#     image(aux, breaks = bks, col=cols, xlab="", ylab="", xaxt=xaxt, yaxt=yaxt,add=add, ...) #yaxt=ifelse(j==1,"s","n")
#     xgrid = seq(from=floor(min(x)), to = ceiling(max(x)), by=1)
#     ygrid = seq(from=floor(min(y)), to = ceiling(max(y)), by=1)
#     abline(lm(y~x), col=2, lwd=2)
#     if(grid)abline(v=xgrid, h=ygrid, col="#999999")
#     legend(legpos, legend=round(cor(x,y), dig=3), bg="#999999")
#   }





#' Swath plots
#' 
#' Plots of data vs. one spatial coordinate, with local average spline curve
#'
#' @param data data to be represented, compositional class, data.frame, or
#' spatial data object (in which case, \code{loc} is a formula!)
#' @param ... further arguments to panel plots
#'
#' @return Nothing, as the function is primarily called to produce a plot
#' @export
#' @importFrom stats loess
#' @importFrom graphics plot par text lines mtext
#'
#' @examples
#' \donttest{
#' data("Windarling")
#' library("sp")
#' compo = compositions::acomp(Windarling[,c("Fe","Al2O3","SiO2", "Mn", "P")])
#' Northing = Windarling$Northing
#' swath(compo, Northing)
#' wind.spdf = sp::SpatialPointsDataFrame(sp::SpatialPoints(Windarling[,6:7]), 
#'      data=compo)
#' swath(wind.spdf, loc=Northing)
#' }
swath <- function(data, ...) UseMethod("swath", data)



#' @describeIn swath swath plot default method
#' @export
#' @method swath default
#' @param loc a numeric vector with the values for one coordinate
#' @param pch symbol to be used for the individual points, defaults to "."
#' @param withLoess either logical (should a loess line be added?) or a list 
#' of arguments to DescTools::lines.loess
#' @param commonScale logical or NA: should all plots share the same vertical 
#' range? FALSE=no, TRUE=yes (default), for
#' compositional data sets, the value NA (=all plots within a row) is also 
#' permitted and is actually the default
#' @param xlab label for the common x axis (defaults to a deparsed version 
#' of loc)
#' @param mfrow distribution of the several plots; it has a good internal default 
#' (not used for compositional 
#' classes)
swath.default <- function(data,  # data (matrix, rmult, aplus, rplus or data.frame)
                           loc,   # X or Y coordinates (a vector)
                           pch=".",
                           withLoess=TRUE, # either logical (should a loess line be added?) or a list of arguments to lines.loess
                           commonScale=TRUE, # logical: should all plots share the same vertical range? FALSE=no, TRUE=yes
                           xlab = deparse(substitute(loc)),
                           ..., # extra arguments to plot
                           mfrow # automatic distribution of figs
){
  # set of swath plots for each possible pwlr, eventually with a loess line 
  if(missing(mfrow)) mfrow = c(floor(sqrt(ncol(data))), ceiling(ncol(data)/floor(sqrt(ncol(data)))))
  if(is(data, "Spatial")) return(swath_SpatialPointsDataFrame(data=data, loc=loc, pch=pch, xlab=xlab, mfrow=mfrow,
                                                    withLoess=withLoess, commonScale=commonScale, ...))
  # preparations
  opar = par()
  opar = par_remove_readonly(opar)
  
  on.exit(par(opar))
  col0 = spectralcolors(10)[10]
  
  # isometric representation
  if(is(data, "data.frame")) data = compositions::rmult(as.matrix(data), V=gsi.getV(data), orig=gsi.orig(data))
  comp = idt(data)
  
  # dimensions
  D = ncol(comp)
  N = length(loc)
  
  # scale the variables
  if(is.na(commonScale)){
    rgY = rep(range(unclass(comp), na.rm=TRUE), D)
    dim(rgY) = c(2,D)
  }else{
    rgY = rep(NA, 2*D) 
    dim(rgY) = c(2,D)
  }
  
  # do the plot
  par(mfrow=mfrow, mar=c(1,1,3,1)/5, oma=c(5,5,5,3))
  for(k in 1:D){
    ij = par()$mfg
    i = ij[1]
    j = ij[2]
    plot(loc, comp[,k], ylim = range(c(comp[,k], rgY[,i]), na.rm=TRUE), pch = pch,
         xaxt = ifelse(i==D,"s","n"), yaxt = ifelse(j==1,"s","n"),...)
    title(main=colnames(comp)[k], line=0.5)
    if(is.logical(withLoess)){
      if(withLoess & requireNamespace("DescTools")){
        a = stats::loess(comp[,k]~loc)
        DescTools::lines.loess(a, col=col0)
      }
    }else if(is(withLoess, "list") & requireNamespace("DescTools")){
      args = withLoess
      if(!("col" %in% names(args))) args$col = col0
      args$x = stats::loess(comp[,k]~loc)
      do.call("DescTools::lines.loess", args=args)
    }
    # add axes on extra sides
    if(i==1) axis(side=3)
    if(j==D) axis(side=4)
  }
  # add labels
  mtext(xlab, side=1, outer=TRUE, line=2.5, cex=1.25)
  # return the old graphical parameters to freeze the plot
  #invisible(opar)
}


#' @describeIn swath  Swath plots for acomp objects
#' @method swath acomp 
#' @export
#' @importFrom compositions acomp
swath.acomp <- function(data,  # composition (rcomp, acomp, ccomp)
           loc,   # X or Y coordinates (a vector)
           pch = ".",
           withLoess = TRUE, # either logical (should a loess line be added?) or a list of arguments to lines.loess
           commonScale = NA, # logical or NA:
           xlab = deparse(substitute(loc)),
           ... # extra arguments to plot
  ){
    # recover byRowCommonScale from commonScale
    byRowCommonScale = ifelse(is.na(commonScale), TRUE, ifelse(commonScale, NA, FALSE))
    # set of swath plots for each possible pwlr, eventually with a loess line 

    # preparations
    opar = par()
    opar = par_remove_readonly(opar)
    
    on.exit(par(opar))
    col0 = spectralcolors(10)[10]
    comp = data
    
    # dimensions
    D = ncol(comp)
    N = length(loc)

    # compute pwlrs
    ij = expand.grid(i=1:D, j=1:D) # pairs of indices
    mat = matrix(0, nrow=D, ncol=D*D)   # matrix of +1 numerator, -1 denominator
    for(k in 1:nrow(ij)){
      mat[ij[k,1],k]=-1
      mat[ij[k,2],k]=1
      if(ij[k,1]==ij[k,2]) mat[ij[k,1],k]=0
    }
    Zpwlr = unclass(as.matrix(cdt(comp))) %*% mat
    
    # scale the variables
    if(is.na(byRowCommonScale)){
      rgY = rep(range(Zpwlr, na.rm=TRUE), 2*D)
      dim(rgY) = c(2,2*D)
    }else{
      if(byRowCommonScale){
        rgY = sapply(1:D, function(i)  range(Zpwlr[,(D*(i-1)+1:D)[-i]], na.rm=TRUE))
      }else{
        rgY = rep(NA, 2*D)
        dim(rgY) = c(2,D)
      }
    }
    
    # do the plot
    par(mfrow=c(D,D), mar=c(1,1,1,1)/5, oma=c(5,5,5,3))
    for(i in 1:D){
      for(j in 1:D){
        k = (i-1)*D+j        
        if(i==j){
          # diagonal plots, only labels
          plot(loc, Zpwlr[,k], type="n", xaxt=ifelse(i==D,"s","n"), yaxt=ifelse(j==1,"s","n") )
          text(mean(range(loc, na.rm=TRUE)), 0, colnames(comp)[i], cex=2)
        }else{
          # off-diagonal plots, swaths  
          plot(loc, Zpwlr[,k], ylim = range(c(Zpwlr[,k], rgY[,i]), na.rm=TRUE), pch=pch,
               xaxt=ifelse(i==D,"s","n"), yaxt=ifelse(j==1,"s","n"), ...)
          if(is.logical(withLoess)){
            if(withLoess){
              a = stats::loess(Zpwlr[,k]~loc)
              DescTools::lines.loess(a, col=col0)
            }
          }else if(is(withLoess, "list")){
            args = withLoess
            if(!("col" %in% names(args))) args$col = col0
            args$x = stats::loess(Zpwlr[,k]~loc)
            do.call("DescTools::lines.swath", args=args)
          }
        }
        # add axes on extra sides
        if(i==1) axis(side=3)
        if(j==D) axis(side=4)
      }
    }
    # add labels
    label1 = c("+ part", "numerator") [1+is.acomp(data)]
    label2 = c("- part", "denominator") [1+is.acomp(data)]
    mtext(label1, side=2, outer=TRUE, line=2.5, cex=1.25)
    mtext(label2, side=3, outer=TRUE, line=2.5, cex=1.25)
    mtext(xlab, side=1, outer=TRUE, line=2.5, cex=1.25)
    # return the old graphical parameters to freeze the plot
    #invisible(opar)
}

#' @describeIn swath Swath plots for ccomp objects
#' @method swath ccomp
#' @export
#' @importFrom compositions ccomp
swath.ccomp <- swath.acomp

#' @describeIn swath Swath plots for rcomp objects
#' @method swath rcomp 
#' @export
#' @importFrom compositions rcomp
swath.rcomp <- swath.acomp




swath_SpatialPointsDataFrame <- function(data, loc, ...){
  swath(data@data, loc=loc, ...)
}


