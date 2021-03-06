% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gstatCompatibility.R,
%   R/compositionsCompatibility.R
\name{fit_lmc}
\alias{fit_lmc}
\alias{fit_lmc.default}
\alias{fit_lmc.logratioVariogramAnisotropy}
\alias{fit_lmc.gstatVariogram}
\alias{fit_lmc.logratioVariogram}
\title{Fit an LMC to an empirical variogram}
\usage{
fit_lmc(v, ...)

\method{fit_lmc}{gstatVariogram}(v, g, model, ...)

\method{fit_lmc}{default}(v, g, model, ...)

\method{fit_lmc}{logratioVariogram}(v, g, model, ...)

\method{fit_lmc}{logratioVariogramAnisotropy}(v, g, model, ...)
}
\arguments{
\item{v}{empirical variogram}

\item{...}{further parameters}

\item{g}{spatial data object, containing the original data}

\item{model}{LMC or variogram model to fit}
}
\value{
Method fit_lmc.gstatVariogram is a wrapper around \code{\link[gstat:fit.lmc]{gstat::fit.lmc()}}, that calls this function
and gives the resulting model its appropriate class (c("variogramModelList", "list")).
Method fit_lmc.default returns the fitted lmc (this function currently uses gstat as a
calculation machine, but this behavior can change in the future)
}
\description{
Fit a linear model of coregionalisation to an empirical variogram
}
\section{Methods (by class)}{
\itemize{
\item \code{gstatVariogram}: wrapper around gstat::fit.lmc method

\item \code{default}: flexible wrapper method for any class for which methods
for \code{\link[=as.gstatVariogram]{as.gstatVariogram()}}, \code{\link[=as.gstat]{as.gstat()}} and \code{\link[=as.variogramModel]{as.variogramModel()}} exist.
In the future there may be direct specialised implementations not depending on
package gstat.

\item \code{logratioVariogram}: method for logratioVariogram wrapping compositions::fit.lmc.
In the future there may be direct specialised implementations,
including anisotropy (not yet possible).

\item \code{logratioVariogramAnisotropy}: method for logratioVariogram with anisotropry
}}

\examples{
data("jura", package = "gstat")
X = jura.pred[,1:2]
Zc = jura.pred[,7:13]
gg = make.gmCompositionalGaussianSpatialModel(Zc, X, V="alr", formula = ~1)
vg = variogram(gg)
md = gstat::vgm(model="Sph", psill=1, nugget=1, range=1.5)
gg = fit_lmc(v=vg, g=gg, model=md)
variogramModelPlot(vg, model=gg)
}
