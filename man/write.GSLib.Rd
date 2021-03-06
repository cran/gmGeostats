% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spSupport.R
\name{write.GSLib}
\alias{write.GSLib}
\title{Write a regionalized data set in GSLIB format}
\usage{
write.GSLib(x, file, header = basename(file))
}
\arguments{
\item{x}{regionalized data set}

\item{file}{filename}

\item{header}{the first line of text for the file, defaults to filename}
}
\value{
The status of closing the file, see \code{\link{close}}
for details, although this is seldom problematic. This function is basically called
for its side-effect of writing a data set in the simplified Geo-EAS format that is
used in GSLIB.
}
\description{
Write a regionalized data set in plain text GSLIB format
}
\examples{
data("jura", package="gstat")
\dontrun{write.GSLib(jura.pred, file="jurapred.txt")}
}
\seealso{
\url{http://www.gslib.com/gslib_help/format.html}
}
