\name{Flamb2}
\alias{Flamb2}

\title{
Lambert projection
}
\description{
  Projects the geo-coordinates into rectangular ones using the Lambert
  projection.
}
\usage{
Flamb2(geoconfig, latrf1 = NA, latrf2 = NA, latref = NA, lngref = NA)
}

\arguments{
  \item{geoconfig}{
    \eqn{n \times 2} matrix, containing geo-coordinates in format
    (\code{lat, -long}, in degrees) of n locations.
  }
  \item{latref, lngref}{
    latitute and -longitude of the reference point (mid-point if not
    provided)
  }
  \item{latrf1, latrf2}{
    range of latitudes used in the projection (default: \code{latref} \eqn{+/- .3 * range(latitudes)})
  }
}

\value{
  A list containing the following named components:
  \item{xy}{\eqn{n \times 2} matrix of Lambert projection coordinates}
  \item{latref, lngref}{}
  \item{latrf1, latrf2}{}
}
