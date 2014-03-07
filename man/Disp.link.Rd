\name{Disp.link}
\alias{Disp.link}

\title{
  Function to plot and link points in a dispersion-distance plot and a
  geographic map.
}
\description{
  User may identify points on the dispersion scatter in order to
  identify line segments on the coordinate plot, **and/or** identify
  individual sites on the coordinate plot in order to highlight the
  corresponding set of points on the dispersion scatter.
}

\usage{
Disp.link(disp.mx, coords, dmap, ddisp, names, device = getOption('device'))
}

\arguments{
  \item{disp.mx}{
    nxn matrix of dispersions (i.e. Var(Z(x)-Z(y))
  }
  \item{coords}{
    nx2 matrix of coordinates
  }
  \item{dmap}{
    device number for existing window to be used for map (optional)
  }
  \item{ddisp}{
    device number for existing window to be used for dispersion plot (optional)
  }
  \item{names}{
    not yet used, but might be used for labelling
  }
  \item{device}{
    a character string giving the name of a function, or the function
    object itself, which when called creates a new graphics device 
  }
}

\value{
  Indices of station pairs selected in the dispersion plot, and indices
  of individual stations selected on the coordinate plot.
}

\note{
  Uses \code{\link{setplot}} to set up coordinates for geographic map
  and \code{\link{Fdist}} to compute distances from coords.
}
