% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gmSpatialModel.R
\name{make.gmMultivariateGaussianSpatialModel}
\alias{make.gmMultivariateGaussianSpatialModel}
\title{Construct a Gaussian gmSpatialModel for regionalized multivariate data}
\usage{
make.gmMultivariateGaussianSpatialModel(
  data,
  coords = attr(data, "coords"),
  model = NULL,
  beta = model$beta,
  formula = model$formula,
  ng = NULL,
  nmax = ng$nmax,
  nmin = ng$nmin,
  omax = ng$omax,
  maxdist = ng$maxdist,
  force = ng$force
)
}
\arguments{
\item{data}{either a data set of any data.frame similar class, or else a \code{\link[sp:SpatialPoints]{sp::SpatialPointsDataFrame()}} containing it}

\item{coords}{the coordinates of the sampling locations, if no SpatialPointsDataFrame was provided}

\item{model}{a variogram model, of any relevant class}

\item{beta}{(see \code{formula}) the coefficients of dependence of the mean of the random field, if these are known; e.g. if \code{formula=~1} constant mean,
and the mean is indeed known, \code{beta} would be a compositional mean; seldom used directly}

\item{formula}{a formula without left-hand-side term, e.g. ~1 or ~Easting+Northing, specifying what do we know of the
dependence of the mean of the random field; this follows the same ideas than in \code{\link[gstat:gstat]{gstat::gstat()}}}

\item{ng}{optional neighborhood information, typically created with \code{\link[=KrigingNeighbourhood]{KrigingNeighbourhood()}}}

\item{nmax}{optional, neighborhood description: maximum number of data points per cokriging system}

\item{nmin}{optional, neighborhood description: minimum number of data points per cokriging system}

\item{omax}{optional, neighborhood description: maximum number of data points per cokriging system per quadrant/octant}

\item{maxdist}{optional, neighborhood description: maximum radius of the search neighborhood}

\item{force}{optional logical, neighborhood description: if not \code{nmin} points are found inside \code{maxdist} radius,
keep searching. This and all preceding arguments for neighborhood definition are borrowed from \code{\link[gstat:gstat]{gstat::gstat()}}}
}
\value{
A "gmSpatialModel" object with all information provided appropriately structured. See \linkS4class{gmSpatialModel}.
}
\description{
Construct a regionalized multivariate data container to be used for Gaussian-based geostatistics: variogram modelling, cokriging and simulation.
}
\examples{
data("jura", package="gstat")
X = jura.pred[,1:2]
Zc = jura.pred[,7:13]
make.gmMultivariateGaussianSpatialModel(data=Zc, coords=X)
}
\seealso{
\code{\link[=SequentialSimulation]{SequentialSimulation()}}, \code{\link[=TurningBands]{TurningBands()}} or \code{\link[=CholeskyDecomposition]{CholeskyDecomposition()}} for specifying the exact
simulation method and its parameters, \link{predict_gmSpatialModel} for running predictions or simulations

Other gmSpatialModel: 
\code{\link{Predict}()},
\code{\link{as.gmSpatialModel}()},
\code{\link{gmSpatialModel-class}},
\code{\link{make.gmCompositionalGaussianSpatialModel}()},
\code{\link{make.gmCompositionalMPSSpatialModel}()}
}
\concept{gmSpatialModel}
