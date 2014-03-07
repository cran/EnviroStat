\name{Falternate3}
\alias{Falternate3}

\title{
  Simultaneous coordinate and variogram estimation
}

\description{
  Perform simultaneous estimation of coords and exponential or gaussian
  variogram by alternating weighted least squares.
}

\usage{
Falternate3(disp, coords, model = 1, a0 = 0.1, t0 = 0.5,
            max.iter = 50, max.fcal = 100, alter.lim = 50,
            tol = 1e-05, prt = 0, dims = 2, lambda = 0,
            ncoords, dev.mon = NULL, verbose = FALSE)
}

\arguments{
  \item{disp}{
    \eqn{n \times n} dispersion matrix, equal to \eqn{2 - 2 (spatial
    correlation matrix)}
  }
  \item{coords}{
    \eqn{n \times 2} coordinate matrix
  }
  \item{model}{
    type of variogram: 1 for exponential or 2 for gaussian
  }
  \item{a0, t0}{
    initial variogram parameter estimates
  }
  \item{max.iter, max.fcal}{
    control parameter for calls to non-linear optimization routines
    (same values used in MDS step and in variogram step) 
  }
  \item{alter.lim}{
    maximum number of iterations of pairs of alternating calls to
    \code{Fmdsfit3} (coordinate estimation) and
    \code{\link{Fvariogfit3}} (variogram fitting)
  }
  \item{tol}{
    convergence criterion for coordinate estimates
  }
  \item{prt}{
    not used
  }
  \item{dims}{
    dimension of multi-dimensional scaling
  }
  \item{lambda}{
    smoothing parameter
  }
  \item{ncoords}{
    \eqn{n \times 2} optional initial coordinates to use if not G-plane
  }
  \item{dev.mon}{
    Function to open the graphics device used for plots monitoring the
    convergence of objective. If \code{NULL}, monitoring plots will not
    be created.
  }
  \item{verbose}{
    if \code{TRUE}, display the results of each iteration to the console
  }
}

\details{
This version permits dimension > 2 for scaling.
In the plotting we'll use a plot symbol proportional to the
third coordinate.
}

\value{
  A list containing the following named components:
  \item{variogfit}{Fitted variogram parameters with new locations}
  \item{ncoords}{\eqn{n\times 2} matrix with coordinates of new locations}
}

\note{
  This version also passes a smoothing parameter to the optimization.
  This parameter probably is not scaled exactly the same as it is
  in sinterp and this has not been investigated yet.  
}

\section{Warning}{
  make sure that coords are scaled reasonably small
  before attempting to compute; otherwise matrix inversion is
  likely not to work in calculation of bending energy matrix.
}
