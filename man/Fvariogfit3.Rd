\name{Fvariogfit3}
\alias{Fvariogfit3}

\title{
  Fit an exponential or gaussian variogram
}

\description{
  Fit an exponential or gaussian variogram with interdistances and
  corresponding dispersion values.
}

\usage{
Fvariogfit3(disp.lt, h.lt, model = 1, a0 = 0.1, t0 = 0.5,
            max.iter = 25, max.fcal = 100, bep = 0, verbose = FALSE)
}

\arguments{
  \item{disp.lt}{
   vector of spatial dispersion values between locations (one value for
   each pair)
  }
  \item{h.lt}{
    vector of corresponding interdistances
  }
  \item{model}{
    type of variogram: 1 for exponential or 2 for gaussian
  }
  \item{a0, t0}{
    initial variogram parameter estimates
  }
  \item{max.iter, max.fcal}{
    control parameters for calls to non-linear optimization routines
  }
  \item{bep}{
    bending energy penalty, described in Sampson-Guttorp method (JASA 1992)
  }
  \item{verbose}{
    if \code{TRUE}, display fitting details to the console
  }
}

\details{
  Exponential variogram is defined as \eqn{a_1 + (2-a_1) (1-exp(-t_0 h))}.
  Gaussian variogram is defined as \eqn{a_1 + (2-a_1) (1-exp(-t_0 h^2))}.
}

\value{
  A list containing the following named components:
  \item{objf}{sum of residual sum of squares, plus bending energy penalty}
  \item{a, t0}{fitted parameters}
  \item{fit}{fitted values of variogram}
}

