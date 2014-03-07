\name{staircase.hyper.est}
\alias{staircase.hyper.est}

\title{
  Estimate hyperparameters of ungauged sites
}

\description{
  This function combines the results from the \code{\link{staircase.EM}} fit and the
  SG method to estimate the hyperparameters associated with the ungauged sites.
}

\usage{
staircase.hyper.est(emfit, covfit, u, p, g, d0 = NULL)
}

\arguments{
  \item{emfit}{
    Output from the \code{staircase.EM} fit
  }
  \item{covfit}{
    The covariance matrix between all locations (with new locations 
    at the beginning). This is an output from the SG fitting
  }
  \item{u}{
    number of new locations
  }
  \item{p}{
    dimension of the multivariate response
  }
  \item{g}{
    number of stations
  }
  \item{d0}{
    (optional) The degrees of freedom for the new locations (ungauged block)
  }
}

\value{
  List with the following elements:
  \item{Delta.0}{
    The degree of freedoms for the new locations. Equal to
    \code{d0}, if given (must be \eqn{> u*p+2}),
    else
	\code{mean(emfit$delta)} if \eqn{> u*p+2};
    \code{u*p+ min(emfit$delta)} otherwise.
  }
  \item{Lambda.0}{
    Conditional variance between new locations given the gauged
    stations
  }
  \item{Xi0.0}{
    the regression slope (Note: \eqn{\tau_{0i} = {\rm kronecker}(\xi_0 ,
      diag(p))})
  }
  \item{H.0}{
    The variance matrix for the rows of \eqn{\tau^{[u]}}
  }
  
  Also all components of the output of the \code{staircase.EM} fit (for
  blocks 1-K).
}

\seealso{
  \code{\link{staircase.EM}}
}
