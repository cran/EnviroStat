\name{Ftransdraw}
\alias{Ftransdraw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Interactively choose smoothing parameter value
}
\description{
  An interactive function showing the fitted variogram and the mapping
  transformation from the geographical space into D-space, allowing the
  user to interactively choose a suitable value for the smoothing
  parameter ("lambda").
}
\usage{
Ftransdraw(disp, Gcrds, MDScrds, gridstr, sta.names, lambda = 0,
           lsq = FALSE, eye, model = 1, a0 = 0.1, t0 = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{disp}{
    \eqn{n \times n} dispersion matrix, equal to \eqn{2-2 (spatial
    correlation matrix)} 
  }
  \item{Gcrds}{
    \eqn{n \ times 2} coordinate matrix (in G-space)
  }
  \item{MDScrds}{
    Coordinates of new locations -- as calculated by
    \code{\link{Falternate3}}
  }
  \item{gridstr}{
    grid coordinates, obtained from \code{\link{Fmgrid}} using
    \code{Gcrds} as input
  }
  \item{sta.names}{
    names of locations; if not provided, locations will be numbered 1 to n
  }
  \item{lambda}{
    initial value of the smoothing parameter
  }
  \item{lsq}{
    logical flag used in Sampson-Guttorp method 
  }
  \item{eye}{
    eye perspective; if not provided, locations will be selected using
    the provided data
  }
  \item{model}{
    type of variogram: 1 for exponential or 2 for gaussian
  }
  \item{a0, t0}{
    initial variogram parameter estimates
  }
}

\value{
  A list containing the following named components:
  \item{Dcrds}{\eqn{n \times 2} cooridnate matrix of new locations in D-space}
  \item{Ddist}{\eqn{n \times n} matrix of interdistances in D-space}
}
\note{
  The graphic device has to be setup before engaging this function, 
  for example using the 'setplot' function:
  \preformatted{par(mfrow=c(1,2))
	temp <- setplot(coords.lamb, ax=TRUE)}

  In window, after entering the function call in R console, click on the
  graphical device to register the curse before proceeding and this
  interactive function will provide instructions for using a different
  smoothing parameter value; ie. enter a new lambda value in R console,
  then click on the graphical device to see results. }

