% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accuracy.R
\name{plot.accuracy}
\alias{plot.accuracy}
\title{Plot method for accuracy curves}
\usage{
\method{plot}{accuracy}(
  x,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xaxs = "i",
  yaxs = "i",
  type = "o",
  col = "red",
  asp = 1,
  xlab = "confidence",
  ylab = "coverage",
  pty = "s",
  main = "accuracy plot",
  colref = col[1],
  ...
)
}
\arguments{
\item{x}{an \code{\link[=accuracy]{accuracy()}} object}

\item{xlim}{graphical parameters, see \code{\link[graphics:plot.default]{graphics::plot.default()}}}

\item{ylim}{graphical parameters, see \code{\link[graphics:plot.default]{graphics::plot.default()}}}

\item{xaxs}{graphical parameters, see \code{\link[graphics:plot.default]{graphics::plot.default()}}}

\item{yaxs}{graphical parameters, see \code{\link[graphics:plot.default]{graphics::plot.default()}}}

\item{type}{graphical parameters, see \code{\link[graphics:plot.default]{graphics::plot.default()}}}

\item{col}{graphical parameters, see \code{\link[graphics:plot.default]{graphics::plot.default()}}}

\item{asp}{graphical parameters, see \code{\link[graphics:plot.default]{graphics::plot.default()}}}

\item{xlab}{graphical parameters, see \code{\link[graphics:plot.default]{graphics::plot.default()}}}

\item{ylab}{graphical parameters, see \code{\link[graphics:plot.default]{graphics::plot.default()}}}

\item{pty}{graphical parameters, see \code{\link[graphics:plot.default]{graphics::plot.default()}}}

\item{main}{graphical parameters, see \code{\link[graphics:plot.default]{graphics::plot.default()}}}

\item{colref}{color for the reference line 1:1}

\item{...}{further graphical parameters to \code{\link[graphics:plot.default]{graphics::plot.default()}}}
}
\value{
Nothing, called to plot the accuracy curve
\eqn{\{(p_i, \pi_i), i=1,2, \ldots, I\}}, where \eqn{\{p_i\}}
are a sequence of nominal confidence of prediction intervals and each \eqn{\pi_i}
is the actual coverage of an interval with nominal confidence \eqn{p_i}.
}
\description{
Plot an accuracy curve out of the outcome of \code{\link[=accuracy]{accuracy()}}.
}
\seealso{
Other accuracy functions: 
\code{\link{accuracy}()},
\code{\link{mean.accuracy}()},
\code{\link{precision}()},
\code{\link{validate}()},
\code{\link{xvErrorMeasures.default}()},
\code{\link{xvErrorMeasures}()}
}
\concept{accuracy functions}
