\name{bgrid}
\alias{bgrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create a bi-orthogonal grid
}
\description{
  Function to create a biorthogonal grid using the thin-plate solution;
  the grid object created can be used with the \code{\link{draw}}
  function to plot the grid.
}
\usage{
bgrid(start, xmat, coef, xlimt, iter.limit = 10, perpc = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{
    coordinates of the starting point to create the grid - usually near
    the center of region
  }
  \item{xmat}{
    \eqn{n \times 2}{n x 2} matrix containing xy coordinates of locations
  }
  \item{coef}{
    coefficients of the thin-plate solution, for instance as returned by
    \code{\link{sinterp}}
  }
  \item{xlimt}{
    vector of 4 elements, representing x- and y-ranges for the grid to
    be created; if not provided, the ranges of \code{xmat} are used.
  }
  \item{iter.limit}{
    limit of iterations
  }
  \item{perpc}{
    parameter to control spacing of the grid
  }
}
\value{
  A list with following named components:
  \item{grid}{coordinates of points in the grid}
  \item{ngrid}{number of points in the grid}
  \item{fldmag}{
    the \dQuote{gradient} index (see Sampson and Guttorp, 1992, for
    details)
  }
  \item{nn}{
    internally set for the maximum number of grid points allowed
  }
}
\references{
  Sampson P.D., and Guttorp, P., "Nonparametric estimation of
  nonstationary spatial covariance structure", in Journal of the
  American Statistical Association, vol 87, pp 108-119, 1992.
}
\seealso{
  \code{\link{draw}}, \code{\link{sinterp}}
}
