#### MAF calculations (2 versions) ----------------
# compute the MAF basis from two matrices of (alr/ilr/clr)-covariance 
gsi.getMAFbasis = function(C0, CT, tol=1e-12){
  svdT = svd(CT)
  A1 = with(svdT, u %*% diag(ifelse(d>tol, 1/sqrt(d), 0)) )
  B = t(A1) %*% C0 %*% A1
  svdB = svd(B)
  o = nrow(svdB$u):1  
  A = A1[,o] %*% svdB$u[o,]
  return(A)
}

  
gsi.computeExplainedVariance <- function(X, Winv){
  rs = sapply(1:ncol(X), function(i) mvar(rmult(outer(X[,i], Winv[i,]))))
  return(rs)
}


#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @param ... generic functionality arguments
#' @export
#' @importFrom stats loadings princomp
Maf = function(x,...) UseMethod("Maf",x)

#' Generalised diagonalisations
#' Calculate several generalized diagonalisations out of a data set and its empirical 
#' variogram
#' @param x a data set, typically of class "data.frame" or of a compositional class 
#' @param vg empirical variogram, of a kind fitting to the data provided
#' @param i a slicer for the variogram, typically this will be one or more indices of 
#' the lag distance to take. %For other options see code{getEmpVariogramSlice}.
#'
#' @return An object extending \code{c("princomp.CLASSOF(x)",}"\code{\link{princomp}}") 
#' with classes "\code{genDiag}", and an extra class argument depending on the 
#' diagonalisation method chosen.
#' 
#' Function \code{Maf} results carry the extra class "\code{maf}", and they correspond
#' to minimum/maximum autocorrelation factors (MAF) as proposed by Switzer and Green 
#' (1984). In this case, the
#' slicer is typically just the index of one lag distance (defaults to i=2). MAF 
#' provides the analytical solution to the joint diagonalisation of two matrices,
#' the covariance of increments provided by the slicing and the conventional covariance
#' matrix (of the idt transformed values, if appropriate). Resulting factors are ordered
#' in decreasing order of spatial continuity, which might produce surprising 
#' scree-plots for those who are used to see a screeplot of a principal component analysis.
#' 
#' Function \code{UWEDGE} (Uniformly Weighted Exhaustive Diagonalization with Gauss 
#' iterations; Tichavsky and Yeredor, 2009) results carry the extra class "\code{uwedge}". 
#' The function 
#' is a wrapper on \code{jointDiag::uwedge} from package \code{jointDiag} (Gouy-Pailler, 2017). 
#' In this case the 
#' slicer is typically just a subset of indices of lag distances to consider 
#' (defaults to the nearest indexes to mininum, maximum and mean lag distances of the 
#' variogram). UWEDGE iteratively seeks for a pair of matrices (a mixing and a
#' demixing matrices) diagonalises the set of matrices \eqn{M_1, M_2, \ldots, M_K} 
#' given, by minimizing the target quantity 
#' \deqn{Q_{uwedge}(A, W) = \sum_{k=1}^K Tr[E_k^t\cdot E_k],}
#' where \eqn{E_k = (W^t\cdot M_k \cdot W- A\cdot \Lambda_k\cdot A^t)} and 
#' \eqn{\Lambda_k = diag(W^t\cdot M_k \cdot W)} is the resulting diagonalized version of 
#' each matrix. Obtained factors are ordered
#' in decreasing order of spatial continuity, which might produce surprising 
#' scree-plots for those who are used to see a screeplot of a principal component analysis.
#' 
#' Function \code{RJD} results carry the extra class "\code{rjd}". The function 
#' is a wrapper on \code{JADE::rjd} (Miettinen, Nordhausen and Taskinen, 2017), 
#' implementing the Rotational joint diagonalisation method (Cardoso and Souloumiac, 1996). 
#' In this case the 
#' slicer is typically just a subset of indices of lag distances to consider 
#' (defaults to the nearest indexes to mininum, maximum and mean lag distances).
#' This algorithm also served for quasi-diagonalising a set of matrices as in UWEDGE,
#' just that in this case the quantity to minimise is the sum of sequares of all off-diagonal
#' elements of \eqn{A^t\cdot M_k\cdot A} for all \eqn{k=1, 2, \ldots K}.
#' 
#' All these functions produce output mimicking \code{\link{princomp}}, i.e. with 
#' elements
#' \describe{
#'   \item{sdev}{contrary to the output in PCA, this contains the square root of the
#'   metric variance of the predictions obtained for each individual factor; this is the
#'   quantity needed for \code{\link{screeplot}} to create plots of explained variance
#'   by factor}
#'   \item{loadings}{matrix of contributions of each (cdt-transformed) original variable to the new factors}
#'   \item{center}{center of the data set (eventually, represented through \code{\link{cdt}}), 
#'   in compositional methods}
#'   \item{scale}{the scalings applied to each original variable}
#'   \item{n.obs}{number of observations}
#'   \item{scores}{the scores of the supplied data on the new factors}
#'   \item{call}{the call to the function (attention: it actually may come much later)}
#' } 
#' and additionally some of the following arguments, in different order
#' \describe{
#'   \item{invLoadings}{matrix of contributions of each factor onto each original variable}
#'   \item{Center}{compositional methods return here the cdt-backtransformed center}
#'   \item{InvLoadings}{compositional methods return here the clr-backtransformed inverse loadings, so that
#'   each column of this matrix can be understood as a composition on itself}
#'   \item{DownInvLoadings}{compositional methods return here the clr-backtransformed "minus inverse loadings", so that
#'   each column of this matrix can be understood as a composition on itself; details in 
#'   \code{\link{princomp.acomp}} }
#'   \item{C1, C2}{Maf returns the two matrices that were diagonalised}
#'   \item{eigenvalues}{Maf returns the generalized eigenvalues of the diagonalisation of C1 and C2}
#'   \item{gof}{UWEDGE returns the values of the goodness of fit criterion across sweeps}
#'   \item{diagonalized}{RJD returns the diagonalized matrices, in an array of (K,D,D)-dimensions, being
#'   D the number of variables in \code{x}}
#'   \item{type}{a string describing which package and which function was used as a workhorse for
#'   the calculation}
#' }
#'  
#' 
#' NOTE: if the arguments you provide to RJD and UWEDGE are not of the appropriate type
#' (i.e. data.frames or equivalent) the default method of these functions will just attempt
#' to call the basis functions JADE:rjd and JointDiag:uwedge respectively. 
#' This will be the case if you provide \code{x} as a "\code{matrix}", or as
#' an "\code{array}". In those cases, the output will NOT be structured as an extension
#' to princomp results; instead they will be native output from those functions.
#' 
#'  
#' @export
#' @method Maf data.frame
#' @aliases genDiag
#' @references Cardoso, J. K. and Souloumiac A. 1996. Jacobi angles for simultaneous
#' diagonalization. SIAM Journal of Matrix Analysis and Applications 17(1), 161-164.
#' 
#' Gouy-Pailler C., 2017. jointDiag: Joint approximate diagonalization of a set of 
#' square matrices. R package version 0.3. https://CRAN.R-project.org/package=jointDiag
#' 
#' Miettinen J., Nordhausen K., and Taskinen, S., 2017. Blind source separation based 
#' on Joint diagonalization in R: The packages JADE and BSSasymp. Journal of Statistical 
#' Software 76(2), 1-31.
#' 
#' Switzer P. and Green A.A., 1984. Min/Max autocorrelation factors for multivariate 
#' spatial imaging, Stanford University, Palo Alto, USA, 14pp.
#' 
#' Tichavsky, P. and Yeredor, A., 2009. Fast approximate joint diagonalization 
#' incorporating weight matrices. IEEE Transactions on Signal Processing 57, 878 – 891.
#' 
#' @family generalised Diagonalisations
#'
#' @examples
#' require("magrittr")
#' require("gstat")
#' require("compositions")
#' data("jura", package="gstat")
#' gs = gstat(id="Cd", formula=log(Cd)~1, locations=~Xloc+Yloc, data=jura.pred) %>% 
#' gstat(id="Co", formula=log(Cd)~1, locations=~Xloc+Yloc, data=jura.pred) %>% 
#' gstat(id="Cr", formula=log(Cr)~1, locations=~Xloc+Yloc, data=jura.pred) %>% 
#' gstat(id="Cu", formula=log(Cu)~1, locations=~Xloc+Yloc, data=jura.pred) %>% 
#' gstat(id="Ni", formula=log(Ni)~1, locations=~Xloc+Yloc, data=jura.pred) %>% 
#' gstat(id="Pb", formula=log(Pb)~1, locations=~Xloc+Yloc, data=jura.pred) %>% 
#' gstat(id="Zn", formula=log(Zn)~1, locations=~Xloc+Yloc, data=jura.pred)
#' vg = variogram(gs)
#' mf = Maf(aplus(jura.pred[, -(1:6)]), vg=vg)
#' mf
#' mf$loadings
#' biplot(mf)
Maf.data.frame <- function(x, vg, i=2,...){
    cl <- match.call()
  C1 = unclass(var(idt(x)))
  mn = colMeans(cdt(x))
  C2 = 0*C1
  ## BEGIN: to change in the future to use gmEVario 
  vg = as.gstatVariogram(vg)
  #dists = unique(vg$dist)
  #vv = vg[abs(dists[i]-vg$dist)<1e-9*dists[i],]
  vv = sapply(split(vg$gamma, vg$id), function(x)x[i])
  vv = vv[order(names(vv))]
  C2[lower.tri(C2,diag=TRUE)]=vv
  dd = diag(diag(C2))
  C2 = unclass(C2 + t(C2) - dd)
  ## END: to change in the future to use gmEVario 
  
  
  dt.centred = scale(x, center = TRUE, scale=FALSE)
  eT = eigen(C1)
  A1 = with(eT, vectors %*% diag(ifelse(values>1e-12, 1/sqrt(values), 0)) ) ## C1^{-0.5}
  B = t(A1) %*% C2 %*% A1 
  eB = eigen(B)
  A2 = eB$vectors
  
  eigenvalues = eB$values
  ord = order(eigenvalues)
  A2 = A2[, ord]
  W = unclass(A1) %*% unclass(A2)
  rownames(W) = colnames(x)
  # do not normalize
  # W = unclass(t(normalize(rmult(t(W))))) 
  W = unclass(W)
  Winv = solve(W)
   colnames(Winv) = colnames(W)
  facscores = unclass(dt.centred) %*% W
  colnames(facscores) = paste("maf", 1:ncol(facscores), sep="")
  sdevs = sqrt(gsi.computeExplainedVariance(facscores, Winv))
  res = list(sdev=sdevs, loadings = W, center=mn, scale=rep(1, ncol(x)),
             n.obs=nrow(x), scores = facscores, C1=C1, C2=C2, eigenvalues = eigenvalues[ord], 
             invLoadings=Winv)
  if(class(x)!="data.frame"){
    res$Center = cdtInv(mn)
    res$InvLoadings = cdtInv(Winv, orig=x)
    res$InvDownLoadings = cdtInv(-Winv, orig=x)
  }
  res$call = cl
  res$type="gmGeostats::MAF"
  class(res)= c("maf","genDiag",paste("princomp", class(x), sep="."), "princomp")
  return(res)
}

#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @method Maf rmult
#' @export
Maf.rmult <- Maf.data.frame

#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @method Maf aplus
#' @export
Maf.aplus <- Maf.data.frame

#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @method Maf rplus
#' @export
Maf.rplus <- Maf.data.frame

#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @method Maf ccomp
#' @export
Maf.ccomp <- Maf.data.frame


#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @method Maf rcomp
#' @export
Maf.rcomp  <- function(x, vg, i=2,...){
  requireNamespace("compositions", quietly=TRUE)
  cl <- match.call()
  C1 = unclass(clrvar2ilr(cov(x)))
  C2 = unclass(clrvar2ilr(variation2clrvar(vg$vg[i,,]))) 
  clr.centred = scale(cdt(x), center = TRUE, scale=FALSE)
  
  eT = eigen(C1)
  A1 = with(eT, vectors %*% diag(ifelse(values>1e-12, 1/sqrt(values), 0)) ) ## C1^{-0.5}
  B = t(A1) %*% C2 %*% A1 
  eB = eigen(B)
  A2 = eB$vectors
  
  eigenvalues = eB$values
  ord = order(eigenvalues)
  A2 = A2[, ord]
  W = unclass(A1) %*% unclass(A2)
  # do not normalize
  #W = unclass(t(normalize(rmult(t(W))))) 
  Wclr = unclass(t(ilr2clr(t(W))))
  rownames(Wclr) = colnames(x)
  Winv = solve(W)
  Wclr.inv = unclass(ilr2clr(Winv))
    colnames(Wclr.inv) = colnames(x)
  facscores = unclass(clr.centred) %*% Wclr
  colnames(facscores) = paste("maf", 1:ncol(facscores), sep="")
  sdevs = sqrt(gsi.computeExplainedVariance(facscores, Wclr.inv))
  res = list(sdev=sdevs, loadings = Wclr, center=mean(cdt(x),...), scale=rep(1, ncol(x)),
             n.obs=nrow(x), scores = facscores, 
             C1=C1, C2=C2, eigenvalues = eigenvalues[ord], invLoadings=Wclr.inv,
             Center = mean(x, ...),
             InvLoadings = cdtInv(Wclr.inv, orig=x),
             InvDownLoadings = cdtInv(-Wclr.inv, orig=x),
             call = cl, type="gmGeostats::MAF")
  kkk = paste("princomp", class(x), sep=".")
  class(res)= c("maf","genDiag", kkk,  "princomp")
  return(res)
} 

#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @method Maf acomp
#' @export
Maf.acomp <- Maf.rcomp



#### UWEDGE calculations ----------

#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @export
UWEDGE = function(x,...){
  if(!requireNamespace("jointDiag", quietly=FALSE)) stop("UWEDGE requires package 'jointDiag' installed")
  UseMethod("UWEDGE",x)
}

#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @method UWEDGE default
#' @export
UWEDGE.default <- function(x, ...) jointDiag::uwedge(M=x, ...)

#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @method UWEDGE acomp
#' @export
UWEDGE.acomp <-  function(x, vg, i=NULL, ...){
  requireNamespace("compositions", quietly=TRUE)
  cl <- match.call()
  vg = as.logratioVariogram(vg) 
  if(is.null(i)){
    i = c(1, floor(c(0.5,1)*dim(vg$vg)[1]))
  }
  V = loadings(princomp(x))
  C1 = clrvar2ilr(cov(x), V = V)
  d = ncol(x)-1
  M = -0.5*apply(vg$vg[c(1,i),,], 1, clrvar2ilr, V=V )
  dim(M)=c(d,d,length(i)+1)
  M[,,1]=C1
  res = jointDiag::uwedge(M)
  
  ord = order(diag(res$B %*% M[,,length(i)] %*% t(res$B)))
  res$B = res$B[ord,]
  W = t(res$B)
  # do not normalize
  # W = unclass(t(normalize(rmult(t(W))))) 
  W = unclass(W)
  ilrx = ilr(x, V=V)
  mns = mean(ilrx)
  ilrxc = ilrx-mns
  facscores = unclass(ilrxc) %*% W
  Wclr = V %*% W
  rownames(Wclr) = colnames(x)
  colnames(Wclr) <- paste("uwedge", 1:d, sep="")
  Wclr.inv = gsiInv(Wclr)
  colnames(Wclr.inv) = colnames(x)
  colnames(facscores) <- paste("uwedge", 1:d, sep="")
  sdevs = sqrt(gsi.computeExplainedVariance(facscores, Wclr.inv))
  out = list(sdev=sdevs, loadings = Wclr, center=mean(cdt(x)), 
             scale=rep(1, ncol(x)),
             n.obs=nrow(x), scores = facscores, 
             gof=res$criter, invLoadings=Wclr.inv,
             Center = mean(x),
             InvLoadings = cdtInv(Wclr.inv, orig=x),
             InvDownLoadings = cdtInv(-Wclr.inv, orig=x),
             call = cl, type="jointDiag::uwedge")
  kkk = paste("princomp", class(x), sep=".")
  class(out)= c("uwedge","genDiag", kkk,  "princomp")
  return(out)
}  

#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @method UWEDGE rcomp
#' @export
UWEDGE.rcomp <- UWEDGE.acomp




### RJD calculation -----------------------

#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @export
RJD = function(x,...){
  if(!requireNamespace("JADE", quietly=FALSE)) stop("RJD requires package 'JADE' installed")
  UseMethod("RJD",x)
}

#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @method RJD default
#' @export
RJD.default <- function(x, ...){
  JADE::rjd(X=x, ...)
} 


#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @method RJD acomp
#' @export
RJD.acomp <- function(x, vg, i=NULL,  ...){
  requireNamespace("compositions", quietly=TRUE)
  cl <- match.call()
  vg = as.logratioVariogram(vg)
  if(is.null(i))
    i = c(1, floor(c(0.5,1)*dim(vg$vg)[1]))
  d = ncol(x)-1
  i = sort(i,decreasing=TRUE)
  V = loadings(princomp(x))
  M = -0.5*apply(vg$vg[i ,,], 1, clrvar2ilr, V=V )
  dim(M)=c(d,d,length(i))
  res = JADE::rjd(M)
  res$B = t(res$V) 
  
  ord = order(diag(res$B %*% M[,,2] %*% t(res$B)))
  res$B = res$B[ord,]
  W = t(res$B)
  # do not normalize
  # W = unclass(t(normalize(rmult(t(W))))) 
  W = unclass(W)
  ilrx = ilr(x, V=V)
  mns = mean(ilrx)
  ilrxc = ilrx-mns
  facscores = unclass(ilrxc) %*% W
  Wclr = V %*% W
  rownames(Wclr) = colnames(x)
  Wclr.inv = t(Wclr)
  colnames(Wclr.inv) = colnames(x)
  sdevs = sqrt(gsi.computeExplainedVariance(facscores, Wclr.inv))
  o = order(sdevs, decreasing = TRUE)
  Wclr = Wclr[,o]
  facscores = facscores[,o]
  Wclr.inv = Wclr.inv[o,]
  colnames(Wclr) <- rownames(Wclr.inv) <- paste("rjd", 1:d, sep="")
  colnames(facscores) <- paste("rjd", 1:d, sep="")
  out = list(sdev=sdevs[o], loadings = Wclr, center=mean(cdt(x)), 
             scale=rep(1, ncol(x)),
             n.obs=nrow(x), scores = facscores, 
             diagonalized = res$D[o,o,], invLoadings=Wclr.inv,
             Center = mean(x),
             InvLoadings = cdtInv(Wclr.inv, orig=x),
             InvDownLoadings = cdtInv(-Wclr.inv, orig=x),
             call = cl, type="JADE::rjd")
  kkk = paste("princomp", class(x), sep=".")
  class(out)= c("rjd","genDiag", kkk,  "princomp")
  return(out)
}  

#' @describeIn Maf.data.frame  Generalised diagonalisations
#' @method RJD rcomp
#' @export
RJD.rcomp <- RJD.acomp


### genDiag methods ----

#' Predict method for generalised diagonalisation objects
#'
#' @param object a generalized diagonalisation object, as obtained from a call to
#' \code{\link{Maf}}, and on the same page, information on the other diagonalisation
#' methods \code{UWEDGE} or \code{RJD}
#' @param newdata a matrix or data.frame of factor scores to convert back to the original
#' scale (default: the scores element from `object`)
#' @param ... not used, kept for generic compatibility
#'
#' @return A data set or compositional object of the nature of the original data
#' used for creating the genDiag object.
#' @include gmAnisotropy.R
#' @export
#' @method predict genDiag
#' @family generalised Diagonalisations
#' @examples
#' data("jura", package="gstat")
#' juracomp = compositions::acomp(jura.pred[, -(1:6)]) 
#' lrvg = logratioVariogram(data=juracomp, loc=jura.pred[,1:2])
#' mf = Maf(juracomp, vg=lrvg)
#' mf
#' biplot(mf)
#' predict(mf) 
#' unclass(predict(mf)) - unclass(juracomp) # predict recovers the original composition
predict.genDiag = function (object, newdata=NULL, ...) {
  requireNamespace("compositions", quietly=TRUE)
  if(is.null(newdata))
    newdata = object$scores
  if("data.frame" %in% class(newdata))
    newdata = as.matrix(newdata)
  Z = newdata %*% unclass(object$invLoadings)
  if("Center" %in% names(object)){
    Z = cdtInv(Z, orig = object$Center) + object$Center
  }else{
    Z = Z + outer(rep(1, nrow(Z)), object$center)
    Z = data.frame(Z)
    colnames(Z) = names(object$center)
  }
  return(Z)
}


#' Colored biplot for gemeralised diagonalisations
#' Colored biplot method for objects of class genDiag
#' @param x a generalized diagonalisation object, as obtained from a call to
#' \code{\link{Maf}} (or to \code{UWEDGE} or \code{RJD}, on the help page of \code{\link{Maf}}).
#' @param choices which factors should be represented? vector of 2 indices; defaults to
#' c(1,2)
#' @param scale deprecated, kept for coherence with \code{link{biplot.princomp}}
#' @param pc.biplot same as the homonimous argument from \code{link{biplot.princomp}}:
#'  boolean, to scale variables down by sqrt(n) and observations up by the same factor. 
#' @param ... further arguments to \code{\link{coloredBiplot}}
#'
#' @return nothing. Function is called exclusively to produce the plot
#' @export
#' @family generalised Diagonalisations
#' @importFrom compositions coloredBiplot
#' @method coloredBiplot genDiag
#' 
#' @references Mueller, Tolosana-Delgado, Grunsky and McKinley (2021) Biplots for 
#' Compositional Data Derived from Generalised Joint Diagonalization Methods. 
#' Applied Computational Geosciences (under review)
#'
#' @examples
#' data("jura", package="gstat")
#' juracomp = compositions::acomp(jura.pred[, -(1:6)]) 
#' lrvg = logratioVariogram(data=juracomp, loc=jura.pred[,1:2])
#' mf = Maf(juracomp, vg=lrvg)
#' mf
#' compositions::coloredBiplot(mf, xlabs.col=as.integer(jura.pred$Rock)+2)
coloredBiplot.genDiag <- function(x, choices = 1:2, scale=0, pc.biplot=FALSE, ...){
  if(length(choices) != 2) stop("length of choices must be 2")
  if(!length(scores <- x$scores))
    stop(gettextf("object '%s' has no scores", deparse(substitute(x))),
         domain = NA)
  lam <- rep(1, along.with=choices)
  if(is.null(n <- x$n.obs)) n <- 1
  lam <- lam * sqrt(n)
  if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
  if(scale != 0) lam <- lam^scale else lam <- 1
  if(pc.biplot) lam <- lam / sqrt(n)
  coloredBiplot.default(t(t(scores[,choices]) / lam),
                         t(x$invLoadings[choices,] * lam), ...)
  invisible()
}
