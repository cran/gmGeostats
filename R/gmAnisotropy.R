

#### Anisotropy ------------------


#' Convert to anisotropy scaling matrix
#'
#' Convert an anisotropy specification to a scaling matrix
#'
#' @param x an matrix to be tagged as anisotropy scaling matrix
#'
#' @return An anisotropy scaling matrix \eqn{A} is such that for any
#' lag vector \eqn{h}, the variogram model turns isotropic in terms
#' of \eqn{u'=h'\cdot A}. This function does not check any special
#' property for this matrix! You should probably be using `anis_GSLIBpar2A()`
#' isntead, and leave `AnisotropyScaling()` for internal uses.
#'
#'
#' @export
#' @family anisotropy
#' @examples
#' ( l = anis_GSLIBpar2A(angles=30, ratios=0.5))
#' ( ll = unclass(l) )
#' AnisotropyScaling(l)
AnisotropyScaling = function(x){
  class(x) = "AnisotropyScaling"
  return(x)
}

#' Convert to anisotropy scaling matrix
#'
#' Convert an anisotropy specification to a scaling matrix
#'
#' @param x an object convertible to an anisotropy scaling matrix; see details
#'
#' @return A matrix \eqn{A} such that for any lag vector \eqn{h}, the variogram model turns
#'  isotropic in terms of \eqn{u'=h'\cdot A}.
#'
#' @details Method `as.AnisotropyScaling.numeric()` expects a vector of two numbers in 2D,
#' or a vector of 5 numbers in 3D. These are in 2D, the azimuth of maximum continuity (in
#' degrees, clockwise from North) and the anisotropy ratio of short/long range. In 3D
#' these are: 1,2) the azimuth and the dip of the direction of maximal continuity; 3) the
#' angle of rotation around the axis of the first direction; 4,5) the anisotropy ratios of
#' the ranges of the second/first and third/first directions of maximal continuity. All angles
#' are given in degrees, all ratios must be smaller or equal to 1. This follows gstat (and hence
#' GSlib) conventions; see gstat::vgm() for details.
#'
#' @family anisotropy
#' @export
#' @aliases as.AnisotropyScaling.numeric
#' @examples
#' ( l = anis_GSLIBpar2A(angles=30, ratios=0.5))
#' ( ll = unclass(l) )
#' as.AnisotropyScaling(ll)
#' @export
as.AnisotropyScaling <- function(x){ UseMethod("as.AnisotropyScaling", x) }

#' @describeIn as.AnisotropyScaling identity method
#' @method as.AnisotropyScaling AnisotropyScaling
#' @export
as.AnisotropyScaling.AnisotropyScaling = function(x){
  x
}


#' @describeIn as.AnisotropyScaling  from a vector of numbers
#' @method as.AnisotropyScaling numeric
#' @export
as.AnisotropyScaling.numeric <- function(x){
  if(length(dim(x))==0){
    if(length(x)==2) return(anis2D_par2A(ratio=x[2], angle=x[1], inv=TRUE))
    if(length(x)==5) return(anis3D_par2A(ratios=x[4:5], angles=x[1:3], inv=TRUE))
    stop("as.AnisotropyScaling.numeric: only works for length=2 1.angle+1.ratio, or lenth=5 3.angles+2.ratios")
  }else{
    warning("as.AnisotropyScaling.numeric: expected vector, but matrix given; attemting an interpretation as AnisotropyScaling or as AnisotropyRangeMatrix")
    if( all(abs(x-t(x))<1e-12) ){
      return(as.AnisotropyScaling.AnisotropyRangeMatrix(x))
    }else{
      return(AnisotropyScaling(x))
    }
  }
}

#' @describeIn as.AnisotropyScaling from an AnisotropicRangeMatrix
#' @method as.AnisotropyScaling AnisotropyRangeMatrix
#' @export
as.AnisotropyScaling.AnisotropyRangeMatrix = function(x){
  AnisotropyScaling(solve(chol(x)))
}



#' Force a matrix to be anisotropy range matrix,
#'
#' Force a matrix M to be considered an anisotropy range matrix, i.e
#' with ranges and orientations,
#' such that \eqn{u = sqrt(h' * M^{-1} * h)} allows to use an isotropic
#' variogram.
#'
#' @param x  matrix simmetric positive definite (i.e. M above)
#' @param checkValidity boolean, should validity be checked?
#'
#' @return the same matrix with a class attribute
#' @family anisotropy
#' @export
AnisotropyRangeMatrix = function(x, checkValidity=TRUE){
  if(checkValidity){
    if(length(dim(x))==3){
      odim = dim(x)
      x = t(sapply(1:(odim[1]), function(i){ AnisotropyRangeMatrix(x[i,,])} ))
      dim(x) = odim
    }else{
      if(nrow(x)!=ncol(x)) stop("AnisotropyRangeMatrix: square matrix needed")
      if(any(abs(x-t(x))>1e-12 )) stop("AnisotropyRangeMatrix: symmetric matrix needed")
      if(any(eigen(x, only.values = TRUE)$values<(-1e-12))) stop("AnisotropyRangeMatrix: positive definite matrix needed")
    }
  }
  class(x) = "AnisotropyRangeMatrix"
  return(x)
}

#' Force a matrix to be anisotropy range matrix,
#'
#' Force a matrix M to be considered an anisotropy range matrix, i.e
#' with ranges and orientations,
#' such that \eqn{u = sqrt(h' * M^{-1} * h)} allows to use an isotropic
#' variogram.
#'
#' @param x  matrix simmetric positive definite (i.e. M above)
#'
#' @return the same anisotropy, specified as `M`
#' @family anisotropy
#' @export
as.AnisotropyRangeMatrix <- function(x){ UseMethod("as.AnisotropyRangeMatrix", x) }


#' @describeIn as.AnisotropyRangeMatrix  Default conversion to anisotropy range matrix
#' @method as.AnisotropyRangeMatrix default
#' @export
as.AnisotropyRangeMatrix.default <- function(x) as.AnisotropyRangeMatrix(as.AnisotropyScaling(x))


#' @describeIn as.AnisotropyRangeMatrix  identity conversion
#' @method as.AnisotropyRangeMatrix AnisotropyRangeMatrix
#' @export
as.AnisotropyRangeMatrix.AnisotropyRangeMatrix <- function(x) x

#' @describeIn AnisotropyRangeMatrix  Convert from AnisotropyScaling
#' @method as.AnisotropyRangeMatrix AnisotropyScaling
#' @export
as.AnisotropyRangeMatrix.AnisotropyScaling <- function(x) AnisotropyRangeMatrix(solve(tcrossprod(x)), checkValidity=FALSE)






#' Produce anisotropy scaling matrix from angle and anisotropy ratios
#'
#' Produce anisotropy matrix (as the transposed of the Cholesky
#' decomposition) from angle and anisotropy ratios
#'
#' @param ratios vector of two values between 0 and 1 giving the anisotropy ratios of
#' medium/largest smallest/largest ranges
#' @param angles as defined in gstat::vgm (and indeed GSLIB). For `anis2D_par2A` 'angle' is the direction of maximum range, i.e. largest spatial continuity, measured clockwise from North
#' @param inv boolean or integer, see `return` for details
#'
#' @return a 3x3 matrix of anisotropy.
#'
#' If `inv=TRUE` (or 1) the output is a matrix `A` such that `norm(h %*% A)`
#' allows to use isotropic variograms, being `h = c(hx, hy, hz)` the lag vectors.
#'
#' If `inv=FALSE` (or 0) the output is a matrix `A` such that `norm(h %*% solve(A))`
#' allows to use isotropic variograms.
#'
#' Other values are meaningless.
#'
#' @family anisotropy
#' @export
  anis_GSLIBpar2A = function(ratios, angles, inv=FALSE){
  if( (length(ratios)==1) & length(angles)==1){
    anis2D_par2A(ratio=ratios, angle=angles, inv=inv)
  }else if( (length(ratios)==2) & length(angles)==3){
    anis3D_par2A(ratios=ratios, angles=angles, inv=inv)
  }else{
    stop("anis.GSLIBpar2A: error, ratios and angles length incompatible for both 2D and 3D")
  }
}





#' @describeIn anis_GSLIBpar2A 2D case
#' @param ratio an anisotropy ratio (min/max range)
#' @param angle direction of maximum range, i.e. largest spatial continuity, measured
#' clockwise from North
#' @export
#' @examples
#' ## ratio=0.5, azimuth 30?? (i.e. direction 60??)
#' A = anis2D_par2A(1, 30)
#' A
#' AAt = A %*% t(A)
#'  #  project the bisector 1:1 (i.e. 45??)
#' (k = c(1,1,0) %*% A)
#' atan2(k[2], k[1]) * 180/pi  # should be 15
#' sqrt(sum(k^2))
#' sqrt( c(1,1,0) %*% AAt %*% c(1,1,0) )
#' A = anis2D_par2A(0.5, 60)
#' rd = 60 * pi/180
#' A
#' A %*% t(A)
#' c(cos(rd), sin(rd),0) %*% A #  should be 1
#' c(-sin(rd), cos(rd),0) %*% A #  should be +/- sqrt(2)
anis2D_par2A = function(ratio, angle, inv=FALSE){
  ## transposed Cholesky decomp of the rotation matrix from angle of maximum
  #     spatial continuity + ratio min/max ranges
  a = angle * pi/180

  # matrix of coordinates of the new basis (in columns)=
  #    =matrix of transformation from new to old coordinate system
  R = matrix(c(sin(a), cos(a), cos(a), -sin(a)), ncol=2)
  # matrix of transformation from old to new coordinate system
  # Rinv = t(R)
  # scaling:
  S = diag(c(1, ratio^((1-2*inv)/2) ))
  # rotation to new coordinates and scaling:
  A = R %*% S # t(Rinv) %*% S # because t(Rinv) = R in orthogonal matrices
  A = rbind(cbind(A,0),c(0,0,1))
  class(A) = "AnisotropyScaling"
  return(A)
}


#' @describeIn anis_GSLIBpar2A 3D case
#' @export
#' @examples
#' c60 = cos(60*pi/180)
#' s60 = sin(60*pi/180)
#' c30 = cos(30*pi/180)
#' s30 = sin(30*pi/180)
#' #  in the new coordinates, 60cwN is (0,1,0)
#' R60p = anis3D_par2A(ratios=c(1,1), angles=c(60,0,0))
#' c(s60, c60, 0) %*% R60p
#' R6030 = anis3D_par2A(ratios=c(1,1), angles=c(60,30,0))
#' # the original X axis is positive on newX and newY, but negative on newZ
#' c(1,0,0) %*% R6030
#' # rotate first direction 60 degrees azimuth, then dip 30degrees upwards
#' c( c60*c30, -s60*c30, s30) %*% R6030
#' (Ranis = anis3D_par2A(ratios=c(0.5,0.25), angles=c(60,30,0)) )
anis3D_par2A = function(ratios, angles, inv=FALSE){
  # A should be a matrix such that for h=[hx hy, hz] lag vectpr
  #   hn = norm( h %*% A[1:ncol(h),]) allows to use the
  #   isotropic correlogram
  ## transposed Cholesky decomp of the rotation matrix from angle of maximum
  #     spatial continuity + ratio min/max ranges
  a = angles * pi/180
  # matrix of coordinates of the new basis (in columns)=
  #    =matrix of transformation from new to old coordinate system
  r = function(a, i=3, cw=TRUE){
    if(i==2) cw=!cw
    m = matrix( c(cos(a), (-1)^cw*sin(a), 0, -(-1)^cw*sin(a), cos(a), 0, 0,0,1), ncol=3)
    ii = 1:3
    jj = order(c(ii[-i],ii[i]))
    return(m[jj,jj])
  }
  R = r(a[1],i=3,cw=T) %*% r(a[2],i=2,cw=T) %*% r(a[3],i=1,cw=F)
  # matrix of transformation from old to new coordinate system
  # Rinv = t(R)
  # scaling:
  S = diag(c(1, ratios^((1-2*inv)/2) ))
  # rotation to new coordinates and scaling:
  A = R %*% S # t(Rinv) %*% S # because t(Rinv) = R in orthogonal matrices
  class(A) = "AnisotropyScaling"
  return(A)
}
