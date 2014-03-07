\name{pred.dist.simul}
\alias{pred.dist.simul}

\title{
  Simulate from the predictive distribution
}

\description{
  Simulate N- replicates from the predictive distribution for
  a given time point (tpt) from 1 to n (length of the data).
}

\usage{
pred.dist.simul(hyperest, tpt, include.obs = T, N = 1)
}

\arguments{
  \item{hyperest}{
    Output from the \code{\link{staircase.hyper.est}} functions, containing estimates of all hyperparameters
  }
  \item{tpt}{
    A specific time point - from 1 to n corresponding to the 
    number of time points from the data set
  }
  \item{include.obs}{
    If TRUE, the observed data for time \code{tpt}, are also returned
  }
  \item{N}{
    Number of replicates
  }
}

\value{
  A matrix with N rows; the number of columns depends on whether the observed data are returned 

  The columns are organized consistent with the observed data 
  (ie. \eqn{u\times p} ungauged blocks, \eqn{g1\times p},
  \eqn{g_2\times p}, \dots)
}

\note{
  This function could be slow if there are missing data at gauged sites 
  correspondind to the selected time point. That is, it is fastest at time points
  corresponding to Block 1 and slower with higher blocks.

}
