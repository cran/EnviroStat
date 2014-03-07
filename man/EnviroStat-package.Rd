\name{EnviroStat-package}
\alias{EnviroStat-package}
\alias{EnviroStat}
\docType{package}

\title{
  Statistical analysis of environmental space-time processes
}
\description{
  EnviroStat provides functions for spatio-temporal modeling of
  environmental processes and designing monitoring networks for
  them based on an approach described in Le and Zidek (2006).
}

\details{
The modeling approach offered by this package has a number of features:

\itemize{
    \item Conditional on knowing the process parameters the
    environmental process is assumed to have (after a suitable
    transformation if necessary) to be a Gaussian random field
    (GRF).

    \item At every spatial location, the process can yield a
    multiplicity of random responses such as air pollutant
    concentrations.

    \item The approach used in the package lies within a Bayesian
    hierarchical modeling framework  However for computational
    expediency empirical shortcuts are made at higher levels of the
    hierarchical setup. Thus for example most hyperparameters are
    fitted using a type II maximum likelihood approach, eliminating
    the need for the the user to specify them. Thus the package can
    handle large fields of monitoring networks, say with 600 or more
    spatial sites.

    \item The approach does not assume a stationary GRF. Instead it
    takes a nonparametric approach where the spatial covariance
    matrix is left completely unspecified and instead endowed with a
    prior distribution with a hypercovariance matrix that can be
    modeled at level two of the hierarchy, making the method quite
    robust against non-stationarity in the random field.

    \item It presents a approach for designing monitoring networks
    based on the well-known warping method of Sampson and Guttorp
    (1992) as developed with Wendy Meiring.

    \item It allows for missing data, providing that these data are
    missing in blocs of time, which after a regional trend is
    fitted, then become exchangeable. For then the blocs of
    residuals can be permuted the get a decreasing or increasing
    staircase pattern in the data matrix something that is required
    in the approach.

    \item It has been empirically assessed in a number of major
    applications and found to yield well calibrated prediction
    intervals. For example, a 95\% interval will cover their
    predictands about 95\% of the time.
}

}
\author{
  Nhu Le, Jim Zidek, Rick White, and Davor Cubranic.

  Fortran code for Sampson-Guttorp estimation authored by Paul D.
  Sampson, Peter Guttorp, Wendy Meiring, and Catherine Hurley.

  Maintainer: Davor Cubranic <cubranic@stat.ubc.ca
}

\references{
  Le, Nhu D. and James V. Zidek. Statistical Analysis of
  Environmental Space-Time Processes. Springer, New York, 2006.
}

\keyword{ package }
%% Optionally other standard keywords, one per line, from file KEYWORDS in
%% the R documentation directory

\seealso{
  See the package vignette for a guided example of complete
  analysis using the package and the manual for details of
  individual functions.
}
