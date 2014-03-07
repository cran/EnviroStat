\name{ozone.NY}
\alias{ozone.NY}
\docType{data}
\title{
New York ozone concentration levels
}
\description{
  The \code{ozone.NY} data frame consists of hourly \eqn{O_3}
  concentration levels (ppb) from nine stations in New York State
  numbered for our purposes from 1 to 9. These data were originally
  downloaded from the EPA's AIRS air quality site (Environmental
  Protection Agency 2013). Furthermore, preliminay analysis
  suggested a sqrt transformation of the original concentrations to
  symmetrize the data distribution.
}

\usage{data(ozone.NY)}
\format{
  Each row of the data set represents a daily record starting at
  April 1, 1995, and ending at September 30, 1995 (183 days).
  Measurements are recorded by station in columns, with four columns
  per station (for hours 8--12); there are nine stations, so there
  are 36 columns of measurements in total.

  The last six stations have no missing observations while stations
  1, 2, 3 have 2616, 2016, and 72 missing hourly observations,
  respectively.

  \describe{
    \item{\code{month}}{month of the measurement, numeric 4--9}
    \item{\code{weekday}}{day of the week, numeric 2--8}
    \item{\code{sqO3.1}--\code{sqO3.4}}{
        square root of the \eqn{O_3} measurement for station 1 for
        hours 8--12 (8am to noon).
    }
    \item{\eqn{\dots}}{}
    \item{\code{sqO3.33}--\code{sqO3.36}}{
        square root of the \eqn{O_3} measurement for station 9 for
        hours 8--12 (8am to noon).
    }
  }
}
\source{
  US Environmental Protection Agency's AIRS air quality site (2013),
  \url{http://www.epa.gov/ttn/airs/airsaqs/}.
}
\keyword{datasets}
