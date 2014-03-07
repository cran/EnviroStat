\name{Fdist}
\alias{Fdist}

\title{
  Interpoint distances
}

\description{
  Function to compute interpoint distances for an \eqn{n \times p}
  coordinate matrix.
}

\usage{
Fdist(crds)
}

\arguments{
  \item{crds}{
    \eqn{n\times p} matrix of \eqn{n} locations, each with \eqn{p}
    coordinates
  }
}

\value{
  \eqn{n\times n} matrix of interdistances; \eqn{(i, j)} element is the
  distance between locations \eqn{i} and \eqn{j}
}
