% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compositionsCompatibility.R, R/variograms.R
\name{logratioVariogram_gmSpatialModel}
\alias{logratioVariogram_gmSpatialModel}
\alias{variogram_gmSpatialModel}
\title{Variogram method for gmSpatialModel objects}
\usage{
logratioVariogram_gmSpatialModel(
  data,
  ...,
  azimuth = 0,
  azimuth.tol = 180/length(azimuth)
)

variogram_gmSpatialModel(object, methodPars = NULL, ...)
}
\arguments{
\item{data}{the data container (see \linkS4class{gmSpatialModel} for details)}

\item{...}{further parameters to \code{\link[gstat:variogram]{gstat::variogram()}}}

\item{azimuth}{which direction, or directions, are desired (in case of directional variogram)}

\item{azimuth.tol}{which tolerance sould be used for directional variograms?}

\item{object}{a gmSpatialModel object containing spatial data.}

\item{methodPars}{(currently ignored)}
}
\value{
Currently the function is just a convenience wrapper on
the variogram calculation functionalities of package "gstat",
and returns objects of class "\code{gstatVariogram}". Check the
help of \code{gstat::variogram} for further information.
In the near future, methods will be created, which will depend on
the properties of the two arguments provided,  \code{object} and
\code{methodPars}.
}
\description{
Compute the empirical variogram of the conditioning data contained in a \linkS4class{gmSpatialModel} object
}
\section{Functions}{
\itemize{
\item \code{logratioVariogram_gmSpatialModel}: logratio variogram method (see \code{\link[=logratioVariogram]{logratioVariogram()}} for details)
}}

