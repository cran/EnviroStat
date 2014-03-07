\name{staircase.EM}
\alias{staircase.EM}

\title{
  Estimate gauged sites hyperparameters
}

\description{
  Estimate \eqn{\cal{H}_g} hyperparameters of the gauged sites using the
  EM algorithm, using the staircase of the missing data to determine the
  default block structure.
}

\usage{
staircase.EM(data, p = 1, block = NULL, covariate = NULL,
             B0 = NULL, init = NULL, a = 2, r = 0.5, verbose = FALSE,
             maxit = 20, tol = 1e-06)
}

\arguments{
  \item{data}{
    data matrix, grouped by blocks each with stations having the same number of 
    missing observations. The blocks are organized in order of decreasing number 
    of missing observations, ie. block 1 has more missing observations
    than block2.

    Default structure:
    \itemize{
      \item Each column represent data from a station; rows are for time
      \item Blocks are decided based on the number of missing observations
    }
  }
  \item{p}{
    number of pollutants measured at each stations. 
    (first p columns of y are for p pollutants from station 1, block 1).
  }
  \item{block}{
    a vector indicating the number of stations in each block - from 1 to K
  }
  \item{covariate}{
    design matrix for covariates created with \code{model.matrix} with \code{as.factor}
  }
  \item{B0}{
    Provided if the hyperparameter \eqn{\beta_0} (B0) is known and not estimated
  }
  \item{init}{
    Initial values for the hyperparameters; output of this function can be used for that
  }
  \item{a}{
    When p=1, the type-II MLE's for delta's are not available. Delta's are assumed to follow
    a gamma distribution with parameters (a,r)
  }
  \item{r}{
    When p=1, the type-II MLE's for delta's are not available. Delta's are assumed to follow
    a gamma distribution with parameters (a,r)
  }
  \item{verbose}{
    flag for writing out the results at each iteration
  }
  \item{maxit}{
    the default maximum number of iterations
  }
  \item{tol}{
    the convergence level.
  }
}

\details{
  The estimated model is as follows:
  \itemize{
    \item \eqn{data \sim  MVN ( z \times \beta , {\rm kronecker}(I, \Sigma) )}
    \item \eqn{\beta \sim MVN (\beta_0 , {\rm kronecker}(F^{-1} , \Sigma ) )}
    \item \eqn{\Sigma \sim GIW (\Theta , \delta )}
  }

  \eqn{\Theta} is a collection of hyperparameters including \eqn{\xi_0,
    \Omega, \Lambda, H^{-1}}.
}

\value{
  A list with following elements:
  \item{Delta}{
    The estimated degrees freedom for each of the blocks (list)
  }
  \item{Omega}{
    The estimated covariance matrix between pollutants
  }
  \item{Lambda}{
    The estimated conditional covariance matrix between stations in each block given
    data at stations in higher blocks (less missing data) - (list)
  }
  \item{Xi0}{
    The estimated slopes of regression between stations in each blocks and those in higher
    blocks (list). Note that \eqn{\tau_{0i} = {\rm kronecker}(\xi_0, diag(p))} - same across stations
    for each pollutants.
  }
  \item{Beta0}{
    Coefficients - assumed to be the same across stations for each pollutant
  }
  \item{Finv}{
    Scale associated with \eqn{\beta_0}
  }
  \item{Hinv}{
    The estimated hyperparameters (list) - inverse of \eqn{H_j}
  }
  \item{Psi}{
    The estimated (marginal) covariance matrix between stations
  }
  \item{block}{
    From input
  }
  \item{data}{
    From input
  }
  \item{covariate}{
    From input
  }
  \item{Lambda.1K}{
    The inverse Bartlett decomposition (eqn 23?)
  }
}

\seealso{
  \code{\link{staircase.hyper.est}}
}
