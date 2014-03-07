\name{ldet.eval}
\alias{ldet.eval}

\title{
  Log determinant of sub-covariances
}

\description{
  Function to calculate the log determinant off all sub-covariance
  matrices of size (k x k) from a covariance matrix.
}

\usage{
ldet.eval(covmat, k, all = FALSE)
}

\arguments{
  \item{covmat}{
    a covariance matrix (ie. non-negative definite, square and symmetric)
  }
  \item{k}{
    dimension of sub-covariance matrices considered
  }
  \item{all}{
    if \code{TRUE}, returns all combinations with corresponding \eqn{log
      |det| }
  }
}
\value{
  \item{coord.sel}{
    The \eqn{k} coordinates having the largest \eqn{log|det|}
  }
  \item{log.det}{
    The \eqn{log|det|} of the submatrix corresponding the \code{coord.sel}
  }
  \item{all.comb}{
    all combinations and their \eqn{log|det|} if \code{all = TRUE}};
  NULL otherwise.

}

\note{
  Setting \code{all = TRUE} may need additionally a large amount of memory and
  so may not work for a large number of combinations!!
}
