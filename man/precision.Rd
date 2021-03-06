% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accuracy.R
\name{precision}
\alias{precision}
\alias{precision.accuracy}
\title{Precision calculations}
\usage{
precision(x, ...)

\method{precision}{accuracy}(x, ...)
}
\arguments{
\item{x}{an object from which precision is to be computed}

\item{...}{generic functionality, not used}
}
\value{
output depends on input and meaning of the function (the term \code{precision}
is highly polysemic)
}
\description{
Precision calculations
}
\section{Methods (by class)}{
\itemize{
\item \code{accuracy}: Compute precision and goodness for accuracy curves, after Deutsch (1997),
using the accuracy curve obtained with \code{\link[=accuracy]{accuracy()}}. This returns a named vector with
two values, one for \code{precision} and one for \code{goodness}.

Mean accuracy, precision and goodness were defined by Deutsch (1997)
for an accuracy curve \eqn{\{(p_i, \pi_i), i=1,2, \ldots, I\}}, where \eqn{\{p_i\}}
are a sequence of nominal confidence of prediction intervals and each \eqn{\pi_i}
is the actual coverage of an interval with nominal confidence \eqn{p_i}.
Out of these values, the mean accuracy (see \code{\link[=mean.accuracy]{mean.accuracy()}}) is computed as
\deqn{ A = \int_{0}^{1} I\{(\pi_i-p_i)>0\} dp,}
where the indicator \eqn{I\{(\pi_i-p_i)>0\}} is 1 if the condition is satisfied and
0 otherwise. Out of it, the area above the 1:1 bisector and under the accuracy
curve is the precision
\eqn{ P = 1-2\int_{0}^{1} (\pi_i-p_i)\cdot I\{(\pi_i-p_i)>0\} dp, }
which only takes into account those points of the accuracy curve where \eqn{\pi_i>p_i}.
To consider the whole curve, goodness can be used
\deqn{G = 1-\int_{0}^{1} (\pi_i-p_i)\cdot (3\cdot I\{(\pi_i-p_i)>0\}-2) dp.}
}}

\seealso{
Other accuracy functions: 
\code{\link{accuracy}()},
\code{\link{mean.accuracy}()},
\code{\link{plot.accuracy}()},
\code{\link{validate}()},
\code{\link{xvErrorMeasures.default}()},
\code{\link{xvErrorMeasures}()}
}
\concept{accuracy functions}
