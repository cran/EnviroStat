\name{Fmgrid}
\alias{Fmgrid}

\title{
  Generate a grid of points
}
\description{
  Function to generate points on a grid. Points are assembled in an
  \eqn{n \times 2} matrix with \code{NA}'s separating series of vertical
  and horizontal lines in the grid.
}
\usage{
Fmgrid(xlim, ylim, xn = 8, xspace, xres, yn = 8, yspace, yres)
}

\arguments{
  \item{xlim}{
    range in x-axis
  }
  \item{ylim}{
    range in y-axis
  }
  \item{xn, yn}{
    the number of vertical and horizontal lines, respectively. This
    parameter is overridden by \code{xspace} and \code{yspace} if
    specified.
  }
  \item{xspace, yspace}{
    the distance between successive vertical and horizontal lines,
    respectively
  }
  \item{xres, yres}{
    the distance between points generated along horizontal and vertical
    lines, respectively; if xres and yres are not specified, then points
    are generated only at the nodes of intersection of the vertical and
    horizontal lines. Note that these nodes appear in duplicate as
    sequences of points are generated first for the vertical lines and
    then for the horizontal lines.
  }
}

\value{
  A list containing the following named components:
  \item{grid}{coordinates of grid points}
  \item{rl.ind}{Indicator of points generated from vertical or horizontal line}
}
