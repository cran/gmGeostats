#### validation specifications --------------

#' Specify a strategy for validation of a spatial model
#'
#' Specify a strategy to validate a spatial model. Currently only 
#' leave-one-out and n-fold cross-validation are available, each specified 
#' by its own function. Leave-one-out takes no parameter.
#'
#' @param nfolds Either, one integer between 2 and the number of hard conditioning data,
#' specifying how many groups do you want to split the data available; or else a factor
#' specifying these groups
#' @param doAll boolean; should each group be used once for validating the model constructed
#' with the remaining groups; else, only the first group will be used for validation, and the other
#' will be used for training.  
#' @param ... ignored
#'
#' @return An object, a list with an appropriate class, controlling the strategy specified. 
#' This can be of class "NfoldCrossValidation" or of class  c("LeaveOneOut", "NfoldCrossValidation").
#' @export
#' @family validation functions
#'
#' @examples
#' NfoldCrossValidation(nfolds=5, doAll=FALSE)
NfoldCrossValidation = function(nfolds=2, doAll=TRUE, ...){
  res = list(strategy="n-fold cross-validation", nfolds=nfolds, doAll=doAll)
  class(res) = "NfoldCrossValidation"
  return(res)
}


#' Specify the leave-one-out strategy for validation of a spatial model
#' 
#' Function to specify the leave-one-out strategy for validation of a spatial model
#' 
#' @return an object of class c("LeaveOneOut", "NfoldCrossValidation") to be used
#' in a call to [validate()]
#' @family validation functions
#' 
#' @export
#' @examples 
#' LeaveOneOut()
LeaveOneOut = function(){
  res = list(strategy="leave-one-out")
  class(res) = c("LeaveOneOut", "NfoldCrossValidation")
  return(res)
}


#### validation methods -------------
#' Validate a spatial model
#' 
#' Validate a spatial model by predicting some values. Typically this will be a validation set,
#' or else some subset of the conditing data. 
#'
#' @param object spatial model object, typically of class [gstat::gstat()] or [gmSpatialModel-class]
#' @param strategy which strategy to follow for the validation? see functions in 'see also' below.
#' @param ... generic parameters, ignored.
#'
#' @return A data frame of predictions (possibly with kriging variances and covariances, or equivalent 
#' uncertainty measures) for each element of the validation set
#' @export
#' @family validation functions 
#' @family accuracy functions
#'
#' @examples
#' data("Windarling")
#' X = Windarling[,c("Easting","Northing")]
#' Z = compositions::acomp(Windarling[,c(9:12,16)])
#' gm = make.gmCompositionalGaussianSpatialModel(data=Z, coords=X)
#' vg = variogram(gm)
#' md = gstat::vgm(range=30, model="Sph", nugget=1, psill=1)
#' gs = fit_lmc(v=vg, g=gm, model=md) 
#' \dontrun{ v1 = validate(gs, strategy=LeaveOneOut()) # quite slow }
#' vs2 = NfoldCrossValidation(nfolds=sample(1:10, nrow(X), replace=TRUE))
#' vs2
#' \dontrun{ v2 = validate(gs, strategy=vs2) # quite slow }
validate <- function(object, strategy, ...) UseMethod("validate", strategy)


#' @describeIn validate Validate a spatial model
#' @method validate LeaveOneOut
#' @export
validate.LeaveOneOut = function(object, strategy, ...){
  if("gstat" %in% class(object)){
    n = nrow(object$data[[1]]$data)
  }else if(is(object, "gmSpatialModel")){
    n = nrow(object@data)
  }else{
    object = try(as.gmSpatialModel(object))
    if(class(object)=="try-error") stop("validate.LeaveOneOut: provided object not interpretable")
    n = nrow(object@data)
  }
  v = validate(object, NfoldCrossValidation(nfolds=n, doAll=TRUE))
  return(v)
}

#' @describeIn validate Validate a spatial model
#' @method validate NfoldCrossValidation
#' @export
validate.NfoldCrossValidation = function(object, strategy, ...){
  if("gstat" %in% class(object)){
    warning("validate: object provided is of class 'gstat', attempting 'gstat.cv(..., remove.all=TRUE, all.residuals=TRUE)'")
    return(gstat::gstat.cv(object, nfold=strategy$nfolds, remove.all = TRUE, all.residuals = TRUE))
  }
  object = try(as.gmSpatialModel(object))
  if(class(object)=="try-error") stop("validate.NfoldCrossValidation: provided object not interpretable")
  # interpret the information about the n-folds provided
  n = strategy$nfolds
  m = nrow(object@data)
  if(length(n)==1){
    nmax = n
    n0 = rep(1:n, each=ceiling(m/n))
    n = sample(x=n0, size=m)      
  }else{
    nmax = max(n)
  }
  if(length(n)!=m) stop("validate.NfoldCrossValidation: provided n-fold info not interpretable (should be either an integer, or a grouping factor)")
  myfun = function(i){
    tk = (n==i)
    object0 = object
    object0@data=object@data[!tk,]
    object0@coords=object@coords[!tk,]
    rs = predict(object0, newdata = object@coords[tk,])
  }
  requireNamespace("foreach", quietly = TRUE)
  if(!strategy$doAll) return(myfun(1))
  i = 1:nmax
  res = foreach(i=i, .combine ="list") %dopar% myfun(i)
  res = unsplit(res, f=n)
  return(res)
}




