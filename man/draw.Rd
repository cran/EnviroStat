\name{draw}
\alias{draw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot the biorthogonal grid
}
\description{
  Function to plot the biorthogonal grid created by the
  \code{\link{bgrid}} function
}
\usage{
draw(data, fs = FALSE, lwidth = c(1, 1), lcolor = c(1, 1),
     cutpts, limits = FALSE, optlist, pts = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    grid object, typically as returned by the \code{\link{bgrid}} function
  }
  \item{fs}{
    indicator of whether non-solid lines can be used
  }
  \item{lwidth}{
    range of line widths to be used.
  }
  \item{lcolor}{
    range for colours
  }
  \item{cutpts}{
    cutpoints; to be created from data if not provided
  }
  \item{limits}{
    indicator if lower and upper limits are added to the cutpts.
  }
  \item{optlist}{
    list containing values of options including \code{lwidth},
    \code{lcolor}, \code{cutpts}, \code{limits}; values are obtained
    from \code{data} if not provided.
  }
  \item{pts}{
    indicator whether points are also plotted.
  }
}
\value{
  Besides drawing the grid, the function also returns \code{optlist}.
}
\references{
  Sampson P.D., and Guttorp, P., "Nonparametric estimation of
  nonstationary spatial covariance structure", in Journal of the
  American Statistical Association, vol 87, pp 108-119, 1992.
}
\note{
  See Sampson+Guttorp (1992) for detail on interpretation of the plotted
  grid, (E.g., solid lines indicate contraction and dashed lines
  indicate expansion.)
}

\seealso{
  \code{\link{bgrid}}
}
