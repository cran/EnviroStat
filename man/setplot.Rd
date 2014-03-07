\name{setplot}
\alias{setplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Set up device for locations
}
\description{
  Function to setup the graphical device to cover the range of locations.
}
\usage{
setplot(xdata, ydata, pretty.call = TRUE, maxdim, axes = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xdata}{
    vector x-axis coordinates of locations, or either a \eqn{n \times 2}
    matrix or a list with named components \code{x} and \code{y} for x-
    and y-axis coordinates, respectively
  }
  \item{ydata}{
    vector y-axis coordinates of locations, not needed if \eqn{n \times
    2} matrix is given for xdata
  }
  \item{pretty.call}{
    if \code{TRUE}, coordinate locations are replaced with a grid
    calculated using \code{\link[base]{pretty}}
  }
  \item{maxdim}{
    plot dimensions, \code{(width, height)}, in inches; if not provided,
    current \code{par("pin")} is used
  }
  \item{axes}{
    logical value indicating whether both axes should be drawn on the
    plot
  }
}
\value{
  A list containing the following named components:
  \item{xlim, ylim}{plot limits for the x- and y-axis, respectively}
  \item{oldpin}{
    old plot dimensions, as returned by
    \code{\link[graphics]{par}('pin')}
  }
  \item{newpin}{new plot dimensions}
}
