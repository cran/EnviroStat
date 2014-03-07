\name{corrfit}
\alias{corrfit}

\title{
  Estimate location correlations
}

\description{
  Function to estimate correlations between all the locations(new+stations)
  using the results of the SG step
}

\usage{
corrfit(crds, Tspline, sg.fit, model = 1)
}

\arguments{
  \item{crds}{
    coordinates of all locations beginning with new locations
  }
  \item{Tspline}{
    the thin-spline fit from the SG-steps
  }
  \item{sg.fit}{
    the mapping resulted from the SG method
  }
  \item{model}{
    variogram model; 1: exponential  2: gaussian
  }
}

\value{
  correlation matrix among the locations
}
