% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gstatCompatibility.R
\name{as.gstat}
\alias{as.gstat}
\alias{as.gstat.default}
\title{Convert a regionalized data container to gstat}
\usage{
as.gstat(object, ...)

\method{as.gstat}{default}(object, ...)
}
\arguments{
\item{object}{regionalized data container}

\item{...}{accessory parameters (currently not used)}
}
\value{
A regionalized data container of class "gstat",
eventually with variogram model included. See \code{\link[gstat:gstat]{gstat::gstat()}} for more info.
}
\description{
Convert a regionalized data container to a "gstat" model object
}
\section{Functions}{
\itemize{
\item \code{as.gstat.default}: default does nothing
}}

\examples{
data("jura", package = "gstat")
X = jura.pred[,1:2]
Zc = jura.pred[,7:13]
gg = make.gmCompositionalGaussianSpatialModel(Zc, X, V="alr", formula = ~1)
as.gstat(gg)
}
