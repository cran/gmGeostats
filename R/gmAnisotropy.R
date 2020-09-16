#### Anisotropy ------------------
#' Check for any anisotropy class
#' 
#' Check that an object contains a valid specification of anisotropy, in any form 
#'
#' @param x object to check
#'
#' @return a logical, TRUE if the object is an anisotropy specification; FALSE otherwise 
#' @export
#'
#' @examples
#' a =  anis2D.par2A(0.5, 30)
#' a
#' is.anisotropySpecification(a)
is.anisotropySpecification = function(x){
  "Anisotropy" %in% class(x)  |  inherits(x, "Anisotropy")
} 


#' Convert to anisotropy scaling matrix
#' 
#' Convert an anisotropy specification to a scaling matrix
#'
#' @param x an object convertible to an anisotropy scaling matrix; see details
#'
#' @return A matrix \eqn{A} such that for any lag vector \eqn{h}, the variogram model turns 
#'  isotropic in terms of \eqn{u=A\cdot h}. 
#'
#' @details Method `as.AnisotropyScaling.numeric()` expects a vector of two numbers in 2D, 
#' or a vector of 5 numbers in 3D. These are in 2D, the azimuth of maximum continuity (in 
#' degrees, clockwise from North) and the anisotropy ratio of short/long range. In 3D 
#' these are: 1,2) the azimuth and the dip of the direction of maximal continuity; 3) the 
#' angle of rotation around the axis of the first direction; 4,5) the anisotropy ratios of 
#' the ranges of the second/first and third/first directions of maximal continuity. All angles
#' are given in degrees, all ratios must be smaller or equal to 1. 
#'    
#' @export
#' @examples
#' as.AnisotropyScaling(c(30,0.5))
as.AnisotropyScaling <- function(x) UseMethod("as.AnisotropyScaling", x)


#' @describeIn as.AnisotropyScaling Convert to anisotropy scaling matrix
#' @method as.AnisotropyScaling AnisotropyScaling
#' @export
as.AnisotropyScaling.AnisotropyScaling = function(x) x


#' @describeIn as.AnisotropyScaling  Convert to anisotropy scaling matrix
#' @method as.AnisotropyScaling numeric
#' @export
as.AnisotropyScaling.numeric = function(x){
  if(length(x)==2) return(anis2D.par2A(ratio=x[2], angle=x[1]))
  if(length(x)==5) stop("as.AnisotropyScaling: 3D from 5-vector values not yet implemented") # return(anis3D.par2A(ratios=x[4:5], angles=x[1:3]))
  stop("as.AnisotropyScaling.numeric: only works for length=2 1.angle+1.ratio, or lenth=5 3.angles+2.ratios")
}


#' Produce anisotropy matrix from angle and anisotropy ratios
#' 
#' Produce anisotropy matrix (in Cholesky decomposition) from angle and anisotropy ratios
#' 
#' @param ratio an anisotropy ratio (min/max range)
#' @param angle direction of maximum range, i.e. largest spatial continuity, measured 
#' counterclockwise from East
#' 
#' @return a 3x3 matrix of anisotropy, in Cholesky form
#'
#' @export
anis2D.par2A = function(ratio, angle){
  ## Cholesky decomp of the rotation matrix from angle of maximum
  #     spatial continuity + ratio min/max ranges
  a = angle * pi/180
  R = matrix(c(cos(a), sin(a), -sin(a), cos(a)), ncol=2)
  S = diag(c(1, sqrt(ratio)))
  A = t(R) %*% S
  A = rbind(cbind(A,0),c(0,0,1))
  class(A) = c("AnisotropyScaling","Anisotropy")
    return(A)
}  ### NOT SURE OF IT
#  anis2D.par2A(1,0)
#  anis2D.par2A(0.5,0)
#  A = anis2D.par2A(0.5, 30)
#    h = rbind(c(0,0), c(2,1), c(1,-2))
#      anish2Dist(h, A = A)
#aa = seq(from=0, to=2*pi, length.out=360)
#A = anis2D.par2A(ratio=1/2, angle=45)
#x = cbind(sin(aa), cos(aa)) %*% A
#plot(x, type="l", asp=1)



