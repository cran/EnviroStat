\name{location.NY}
\alias{location.NY}
\docType{data}
\title{
  Locations of New York ozone measuring stations.
}
\description{
  The \code{location.NY} data frame consists of geographical
  coordinates of nine stations in New York State
  numbered for our purposes from 1 to 9.
}

\usage{data(location.NY)}
\format{
  A data frame with 9 observations on the following 2 variables.
  \describe{
    \item{\code{lat}}{Station latitude}
    \item{\code{long}}{Station longitude}
  }
}
\seealso{
  This dataset accompanies the ozone measurements in
  \link{ozone.NY}.
}
\examples{
data(location.NY)
plot(location.NY[,1], location.NY[,2])
}
\keyword{datasets}
