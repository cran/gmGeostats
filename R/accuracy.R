##### accuracy and precision -----------

#' Compute accuracy and precision
#' 
#' Computes goodness-of-fit measures (accuracy, precision and joint goodness) adapted or extended from the 
#' definition of Deutsch (1997).
#'
#' @param object data container for the predictions (plus cokriging error variances/covariance) or simulations (and eventually for the true values in univariate problems)
#' @param ... generic functionality, currently ignored
#'
#' @return If `outMahalanobis=TRUE` (the primary use), this function returns a two-column dataset of class
#' c("accuracy", "data.frame"), which first column gives the nominal probability cutoffs used, and the second column
#' the actual coverage of the intervals of each of these probabilities. If `outMahalanobis=FALSE`, the output
#' is a vector (for prediction) or matrix (for simulation) of Mahalanobis error norms.
#' 
#' @details For method "kriging", `object` must contain columns with names including the string "pred" for predictions
#' and "var" for the kriging variance; the observed values can also be included as an extra column with name "observed",
#' or else additionally provided in argument `observed`. For method "cokriging", the columns of `object` must contain
#' predictions, cokriging variances and cokriging covariances in columns including the strings "pred", "var" resp. "cov",
#' and observed values can only be provided via `observed` argument. Note that these are the natural formats when
#' using [gstat::predict.gstat()] and other (co)kriging functions of that package.
#' 
#' For univariate and multivariate cokriging results (methods "kriging" and "cokriging"), the coverage values are computed based on the
#' Mahalanobis square error, the (square) distance between prediction and true value, using as the positive definite bilinear form 
#' of the distance the variance-covariance cokriging matrix. The rationale is that, under the assumption
#' that the random field is Gaussian, the distribution of this Mahalanobis square error should
#' follow a \eqn{\chi^2(\nu)} with degrees of freedom \eqn{\nu} equal to the number of variables. Having this
#' reference distribution allows us to compute confidence intervals for that Mahalanobis square error, and then 
#' count how many of the actually observed errors are included on each one of the intervals (the *coverage*).
#' For a perfect adjustment to the distribution, the plot of coverage vs. nominal confidence (see [plot.accuracy])
#' should fall on the \eqn{y=x} line. NOTE: the original definition of Deutsch (1997) for univariate case
#' did not make use of the \eqn{\chi^2(1)} distribution, but instead derived the desired intervals (symmetric!) 
#' from the standard normal distribution appearing by normalizing the residual with the kriging variance; the result is the
#' same. 
#' 
#' For method "simulation" and object `object` is a data.frame, the variable names containing the realisations must 
#' contain the string "sim", and `observed` must be a vector with as many elements as rows has `object`. If 
#' `object` is a [DataFrameStack()], then it is assumed that the stacking dimension is running through the realisations;
#' the true values must still be given in `observed`.
#' In both cases, the method is based on ranks:
#' with them we can calculate which is the frequency of simulations being more extreme than the observed value.
#' This calculation is done considering bilateral intervals around the median of (realisations, observed value) 
#' for each location separately. 
#' 
#' Method "mahalanobis" ("Mahalanobis" also works) is the analogous for multivariate simulations. It
#' only works for `object` of class [DataFrameStack()], and requires the stacking dimension to run through
#' the realisations and the other two dimensions to coincide with the dimensions of `observed`, i.e.
#' giving locations by rows and variables by columns. In this case, a covariance matrix will be computed
#' and this will be used to compute the Mahalanobis square error defined above in method "cokriging": 
#' this Mahalanobis square error will be computed for each simulation and for the true value. 
#' The simulated Mahalanobis square errors will then be used to generate the reference distribution 
#' with which to derive confidence intervals. 
#' 
#' Finally, highly experimental "flow" method requires the input to be in the same shape as method 
#' "mahalanobis". The method is mostly the same, just that before the Mahalanobis square errors
#' are computed a location-wise flow anamorphosis ([ana()]) is applied to transform the realisations (including
#' the true value as one of them) to joint normality. The rest of the calculations are done as if with
#' method "mahalanobis".
#' 
#' @export
#' @family accuracy functions
#' @references Mueller, Selia and Tolosana-Delgado (2023) Multivariate cross-validation
#' and measures of accuracy and precision. 
#' Mathematical Geosciences (under review).
#' 
#' 
accuracy <- function(object,...) UseMethod("accuracy", object)



#' @describeIn accuracy Compute accuracy and precision
#' @method accuracy data.frame
#' @param observed either a vector- or matrix-like object of the true values
#' @param prob sequence of cutoff probabilities to use for the calculations
#' @param method which method was used for generating the predictions/simulations? 
#' one of c("kriging", "cokriging", "simulation") for `object` of class "data.frame", or of
#' c("simulation", "mahalanobis", "flow") for `object` of class [DataFrameStack()].
#' @param outMahalanobis if TRUE, do not do the final accuracy calculations and return the Mahalanobis
#' norms of the residuals; if FALSE do the accuracy calculations
#' @param ivar if `method`="kriging" or "cokriging" you can also specify here one single variable name 
#' to consider for univariate accuracy; this variable name must exist both in `object` 
#' (including "pred" and "var" prefixes or suffixes in the column names) and in `observed`;
#' this might require renaming the columns of these files!
#' @export
accuracy.data.frame <- function(object, observed=object$observed, 
                                prob = seq(from=0, to=1, by=0.05),
                                method="kriging", outMahalanobis=FALSE, 
                                ivar, ...){
  x = object # after v 0.11.9003, 'accuracy' first argument is renamed to 'object' for compatibility with tidymodels 
  methods = c("kriging", "cokriging", "simulation")
  mm = methods[pmatch(method, methods)]
  if(!missing(ivar)){
    iTrue = grep(ivar, colnames(observed))
    iPred = intersect(grep(ivar, colnames(x)), grep("pred", colnames(x)))
    iVar = intersect(grep(ivar, colnames(x)), grep("var", colnames(x)))
    if(any(sapply(list(iTrue, iPred, iVar), length)!=1))
      stop("accuracy: univariate case by specifying an `ivar` requires the named variable to occur ONCE in `observed` and once in `object`, here prefixed or suffixed with `pred` and `var` to identify kriging predictions and variances")
    mm = "kriging"
    observed = observed[,iTrue]
    x = x[, c(iPred, iVar)]
  } 
  if(mm=="kriging"){
    mynorms = xvErrorMeasures(x, observed=observed, output = "Mahalanobis", univariate=TRUE)
    D = 1
  }else if(mm=="cokriging"){
    mynorms = xvErrorMeasures(x, observed=observed, output = "Mahalanobis", univariate=FALSE)
    D = ncol(observed)
  }else if(mm=="simulation"){
    oneAcc.sim = function(sims, true){
      rks = rank(c(true,sims))
      2*abs(rks[1]/(1+length(sims))-0.5)
    }
    ids = grep("sim", colnames(x))
    if(outMahalanobis)
      warning("accuracy: outMahalanobis=TRUE not valid with method='simulation'")
    if(length(ids)!=ncol(x)){
      warning("accuracy: x should include only simulations, without spatial coordinates; attempting a patch")
      x = x[,ids]
    }
    if(nrow(x)!=length(observed))
      stop("accuracy: dimensions of x and observed do not match")
    sims = as.matrix(x)
    res = sapply(1:length(observed), function(i) oneAcc.sim(sims[i,], observed[i]))
    aa = outer(res, prob, "<=")
    a = colMeans(aa)
    erg = data.frame(p=prob, accuracy=a)
    class(erg) = c("accuracy", "data.frame")
    return(erg)
  }else
    stop('accuracy: method must be one of c("kriging", "cokriging", "simulation")')
  if(outMahalanobis)
    return(mynorms)
  # cases for kriging and cokriging
  qqq = stats::qchisq(prob,df=D)   # TODO: this could be Hotelling's T^2 distributed? 
  aa = outer(mynorms,qqq,"<=")
  a = colMeans(aa)
  erg = data.frame(p=prob, accuracy=a)
  class(erg) = c("accuracy", "data.frame")
  return(erg)
}



#' @describeIn accuracy Compute accuracy and precision
#' @param ivars in multivariate cases, a vector of names of the variables to analyse (or one single variable name)
#' @method accuracy DataFrameStack
#' @export
accuracy.DataFrameStack <- function(object, observed, 
                                    ivars = intersect(colnames(observed), dimnames(object)[[noStackDim(object)]]),
                                    prob = seq(from=0, to=1, by=0.05),
                                    method = ifelse(length(ivars)==1, "simulation", "Mahalanobis"),
                                    outMahalanobis=FALSE, ...){
  x = object  # after v 0.11.9003, 'accuracy' first argument is renamed to 'object' for compatibility with tidymodels 
  methods = c("simulation", "Mahalanobis","mahalanobis", "flow")
  mm = methods[pmatch(method, methods)]
  oneAcc.sim = function(sims, true){
    rks = rank(c(true,as.matrix(sims)))
    2*abs(rks[1]/(1+length(sims))-0.5)
  }
  oneAcc.assim = function(sims, true){
    rks = rank(c(true,as.matrix(sims)))
    rks[1]/(1+length(sims))
  }
  
  if(mm=="simulation"){
    if(outMahalanobis)
      warning("accuracy: outMahalanobis=TRUE not valid with method='simulation'")
    if(nrow(x)!=length(observed))
      if(nrow(x)!=nrow(observed))
        stop("accuracy: dimensions of `object` and `observed` do not match")
    sims = as.matrix(gmApply(x, FUN=function(xx)xx[,ivars, drop=FALSE]))
    res = sapply(1:nrow(sims), function(i) oneAcc.sim(sims[i,], observed[i, ivars, drop=FALSE]))
    aa = outer(res, prob, "<=")
    a = colMeans(aa)
    erg = data.frame(p=prob, accuracy=a)
    class(erg) = c("accuracy", "data.frame")
    return(erg)
  }
  sim_array = as.array(x)
  permidx = c(ifelse(is.numeric(stackDim(x)),1,names(dimnames(x))[1]),stackDim(x),noStackDim(x)) 
  sim_array = aperm(sim_array, perm = permidx)
  sim_array = sim_array[,,ivars, drop=FALSE]
  observed = as.matrix(unclass(observed)[,ivars, drop=FALSE])
  N = dim(sim_array)[1]
  S = dim(sim_array)[stackDim(x)]
  D = dim(sim_array)[noStackDim(x)]
  if(method=="flow"){
    sim_trafos = gmApply(sim_array, 1, ana)
    xx = lapply(1:N, function(k) sim_trafos[[k]](x=rbind(observed[k,], sim_array[k,,])) )
    observed = t(sapply(xx, function(x)x[1,]))
    sim_array = as.array(sapply(xx, function(x)x[-1,]))
    dim(sim_array)=c(S,D,N)
    sim_array = aperm(sim_array, c(3,1,2))
    method = "Mahalanobis"
  }
  if(method=="Mahalanobis" | method=="mahalanobis"){
    invcovs = lapply(1:N,function(k) gsi.powM(cov(sim_array[k,,]),-1)) 
    means = lapply(1:N,function(k) observed[k,]) 
    if(!outMahalanobis)
      means = lapply(1:N,function(k) colMeans(sim_array[k,,])) 
    myfun = function(k){
      xx = sweep(rbind(observed[k,], sim_array[k,,]),2,means[[k]])
      nmsq = sapply(1+0:S, function(j) xx[j,] %*% invcovs[[k]] %*% xx[j,] )
      return(nmsq)
    }
    mynorms = sapply(1:N, myfun)
    if(outMahalanobis)
      return(t(mynorms)[,-1])
    res = sapply(1:N, function(i) oneAcc.assim(mynorms[-1,i], mynorms[1,i]))
  }else stop("accuracy: method given not existing")
  # stuff to use or remove later:

  # cases for kriging and cokriging
  aa = outer(res, prob, "<=")
  a = colMeans(aa)
  erg = data.frame(p=prob, accuracy=a)
  class(erg) = c("accuracy", "data.frame")
  return(erg)
}


#' Mean accuracy
#' 
#' Mean method for accuracy curves, after Deutsch (1997) 
#'
#' @param x accuracy curve obtained with [accuracy()]
#' @param ... generic functionality, not used
#'
#' @return The value of the mean accuracy after Deutsch (1997); details
#' can be found in \code{\link[=precision.accuracy]{precision()}}.
#' @export
#' @family accuracy functions
mean.accuracy = function(x, ...){
  aux = x$accuracy - x$p
  n = nrow(x)
  n/(n-1)*mean(ifelse(aux>0,1,0),...)
}



#' Precision calculations
#' 
#' Precision calculations 
#'
#' @param x an object from which precision is to be computed
#' @param ... generic functionality, not used
#' @return output depends on input and meaning of the function (the term `precision` 
#' is highly polysemic)
#' @export
precision <- function(x,...) UseMethod("precision",x)


#' @describeIn precision Compute precision and goodness for accuracy curves, after Deutsch (1997),
#' using the accuracy curve obtained with [accuracy()]. This returns a named vector with 
#' two values, one for `precision` and one for `goodness`.
#' 
#' Mean accuracy, precision and goodness were defined by Deutsch (1997) 
#' for an accuracy curve \eqn{\{(p_i, \pi_i), i=1,2, \ldots, I\}}, where \eqn{\{p_i\}}
#' are a sequence of nominal confidence of prediction intervals and each \eqn{\pi_i}
#' is the actual coverage of an interval with nominal confidence \eqn{p_i}.
#' Out of these values, the mean accuracy (see [mean.accuracy()]) is computed as
#' \deqn{ A = \int_{0}^{1} I\{(\pi_i-p_i)>0\} dp,} 
#' where the indicator \eqn{I\{(\pi_i-p_i)>0\}} is 1 if the condition is satisfied and 
#' 0 otherwise. Out of it, the area above the 1:1 bisector and under the accuracy 
#' curve is the precision
#' \eqn{ P = 1-2\int_{0}^{1} (\pi_i-p_i)\cdot I\{(\pi_i-p_i)>0\} dp, }
#' which only takes into account those points of the accuracy curve where \eqn{\pi_i>p_i}.
#' To consider the whole curve, goodness can be used
#' \deqn{G = 1-\int_{0}^{1} (\pi_i-p_i)\cdot (3\cdot I\{(\pi_i-p_i)>0\}-2) dp.}
#' 
#' @family accuracy functions
#' @method precision accuracy
#' @export
precision.accuracy <- function(x, ...){
  aux = x$accuracy - x$p
  a = ifelse(aux>0,1,0)
  erg = c(precision=1-2*(a %*% aux)/(nrow(x)-1),
      goodness = 1-((3*a-2) %*% aux)/(nrow(x)-1))
  return(erg)
}



#' Plot method for accuracy curves
#' 
#' Plot an accuracy curve out of the outcome of [accuracy()].
#'
#' @param x an [accuracy()] object
#' @param xlim graphical parameters, see [graphics::plot.default()]
#' @param ylim graphical parameters, see [graphics::plot.default()] 
#' @param xaxs graphical parameters, see [graphics::plot.default()] 
#' @param yaxs  graphical parameters, see [graphics::plot.default()]
#' @param type  graphical parameters, see [graphics::plot.default()]
#' @param col  graphical parameters, see [graphics::plot.default()]
#' @param asp  graphical parameters, see [graphics::plot.default()]
#' @param xlab  graphical parameters, see [graphics::plot.default()]
#' @param ylab  graphical parameters, see [graphics::plot.default()]
#' @param pty  graphical parameters, see [graphics::plot.default()]
#' @param main  graphical parameters, see [graphics::plot.default()]
#' @param colref  color for the reference line 1:1
#' @param ... further graphical parameters to [graphics::plot.default()]
#'
#' @return Nothing, called to plot the accuracy curve 
#' \eqn{\{(p_i, \pi_i), i=1,2, \ldots, I\}}, where \eqn{\{p_i\}}
#' are a sequence of nominal confidence of prediction intervals and each \eqn{\pi_i}
#' is the actual coverage of an interval with nominal confidence \eqn{p_i}. 
#' @export
#' @method plot accuracy
#' 
#' @family accuracy functions
plot.accuracy <- function(x, xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", type="o", col="red", asp=1,
                          xlab="confidence", ylab="coverage", pty="s", main="accuracy plot",
                          colref=col[1], ...){
  par(pty=pty)
  plot(x$p, x$accuracy, xaxs=xaxs, yaxs=yaxs, xlim=xlim, ylim=ylim, type=type, col=col, asp=asp,
       xlab=xlab, ylab=ylab, main=main, ...)
  abline(a=0, b=1, col=colref)
}


### error measures ----------


#' @describeIn xvErrorMeasures.data.frame Cross-validation errror measures
#' @param ... extra arguments for generic functionality
#' @export
xvErrorMeasures <- function(x,...) UseMethod("xvErrorMeasures", x)



#' Cross-validation errror measures
#' 
#' Compute one or more error measures from cross-validation output
#'
#' @param x a vector containing the predicted values
#' @param krigVar a vector containing the kriging variances
#' @param observed a vector containing the true values
#' @param output  which output do you want? a vector of one or several of  c("ME","MSE","MSDR","Mahalanobis")
#' @param ... extra arguments for generic functionality
#'
#' 
#' @export
#' @method xvErrorMeasures default
#' @family accuracy functions
xvErrorMeasures.default <- function(x, krigVar, observed, output="MSDR1", ...){
  if(length(output)>1)
    return(sapply(output, function(o) xvErrorMeasures(x, krigVar, observed, output=o)))
  
  outputs = c("ME","MSE","MSDR","MSDR1","MSDR2","Mahalanobis")
  output = outputs[pmatch(output, outputs, duplicates.ok = TRUE)]
  
  resids = x - observed
  mynorms = resids^2/krigVar
  if(output=="ME") return(mean(resids, na.rm=TRUE))
  if(output=="MSE") return(mean(resids^2, na.rm=TRUE))
  if(output=="MSDR") return(mean(mynorms, na.rm=TRUE))
  if(output=="Mahalanobis") return(mynorms)
}

#' Cross-validation errror measures
#' 
#' Compute one or more error measures from cross-validation output
#'
#' @param x a dataset of predictions (if `x` is of class "data.frame") or simulations
#'  (if `x` is of class "DataFrameStack")
#' @param observed a vector (if univariate) or a matrix/dataset of true values
#' @param output which output do you want? a vector of one or several of  c("ME","MSE","MSDR","MSDR1","MSDR2","Mahalanobis")
#' @param univariate logical control, typically you should not touch it
#' @param ... extra arguments for generic functionality
#'
#' @return If just some of c("ME","MSE","MSDR","MSDR1","MSDR2") are requested, the output is a named
#' vector with the desired quantities. If only "Mahalanobis" is requested, the output is a vector 
#' of Mahalanobis square errors. If you mix up things and ask for "Mahalanobis" and some of
#' the quantities mentioned above, the result will be a named list with the requested quantities. 
#' (NOTE: some options are not available for `x` a "DataFrameStack")
#' 
#' @details "ME" stands for *mean error* (average of the differences between true values and predicted values), 
#' "MSE" stands for *mean square error* (average of the square differences between true values and predicted values), 
#' and "MSDR" for *mean squared deviation ratio* (average of the square between true values and predicted values
#' each normalized by its kriging variance). These quantities are classically used in evaluating 
#' output results of validation excercices of one single variable.
#' For multivariate cases, "ME" (a vector) and "MSE" (a scalar) work as well, 
#' while two different definitions of a multivariate 
#' mean squared deviation ratio can be given:
#' * "MSDR1" is the average Mahalanobis square error (see [accuracy()] for explanations)
#' * "MSDR2" is the average univariate "MSDR" over all variables.
#' 
#' @export
#' @method xvErrorMeasures data.frame
#' @family accuracy functions
xvErrorMeasures.data.frame = function(x, observed=x$observed, output="MSDR1",
                           univariate=length(dim(observed))==0, ...){
  if(length(output)>1)
    return(sapply(output, function(o) xvErrorMeasures(x=x, observed=observed, univariate=univariate, output=o)))
  
  outputs = c("ME","MSE","MSDR","MSDR1","MSDR2","Mahalanobis")
  output = outputs[pmatch(output, outputs, duplicates.ok = TRUE)]
  
  if(univariate){
    prd = grep("pred", colnames(x))
    vr = grep("var", colnames(x))
    if(any(prd!=1, vr!=2, is.null(observed)))
      warning("xvErrorMeasures: method='univariate' expects x to be the output of a gstat.cv call")
    resids = x[,prd] - observed
    mynorms = resids^2/x[,vr]
    if(output=="ME") return(mean(resids, na.rm=TRUE))
    if(output=="MSE") return(mean(resids^2, na.rm=TRUE))
    if(output=="MSDR" | output=="MSDR1" |output=="MSDR2") return(mean(mynorms, na.rm=TRUE))
    if(output=="Mahalanobis") return(mynorms)
  }else{
    xreord = gsi.gstatCokriging2rmult(x)
      covs = attr(xreord, "krigVar")
    preds = sapply(unclass(xreord), cbind)
    obs = as.matrix(observed)
    resids = preds-obs
    if(output=="ME") return(colMeans(resids, na.rm=TRUE))
    if(output=="MSE") return(mean(rowSums(resids^2, na.rm=T), na.rm=TRUE))
    if(output=="MSDR2"){
      myvar = function(j) covs[,j,j]
      erg = resids^2/sapply(1:ncol(resids), myvar)
      return(mean(colMeans(erg, na.rm=T), na.rm=TRUE))
    }
    myfun = function(i){
      resids[i,] %*% solve(covs[i,,], resids[i,])
    }
    mynorms = sapply(1:nrow(resids), myfun)
    if(output=="MSDR1") return(mean(mynorms, na.rm=TRUE))
    if(output=="Mahalanobis") return(mynorms)
  }
}

#' @describeIn xvErrorMeasures.data.frame Cross-validation errror measures
#' @method xvErrorMeasures DataFrameStack
#' @export
xvErrorMeasures.DataFrameStack = function(x, observed, output="ME",
                                      univariate=length(dim(observed))==0,
                                      ...){
  if(length(output)>1)
    return(sapply(output, function(o) xvErrorMeasures(x=x, observed=observed, univariate=univariate, output=o)))
  
  outputs = c("ME","MSE")
  output = outputs[pmatch(output, outputs, duplicates.ok = TRUE)]
  dims = noStackDim(x)
  dims = c(ifelse(is.numeric(dims), 1, "loc"), dims)
  mn = gmApply(x, MARGIN=dims, FUN=mean)
  # class(mn) = "data.frame"
  # xvErrorMeasures(mn, observed, output, univariate=TRUE)
  
  resids = mn-observed
  myfun = function(output){
    if(output=="ME") return(colMeans(resids, na.rm=TRUE))
    if(output=="MSE") return(mean(rowSums(resids^2, na.rm=T), na.rm=TRUE))
  }
  
  return(sapply(outputs, myfun))
}
