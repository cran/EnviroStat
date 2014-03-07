\name{seval}
\alias{seval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Value estimation using thin-plate spline
}
\description{
  Function to estimate value at a location using thin-plate spline
}
\usage{
seval(x, tpsp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    \eqn{n \times 2}{n x 2} matrix, containing xy coordinates of
    locations
  }
  \item{tpsp}{
    thin-plate solution, typically returned by \code{\link{sinterp}}
  }
}
\value{
  List containing the following named components:
  \item{x}{location coordinates, copied from the \code{x}}
  \item{y}{estimated values at location \code{x}}
}
