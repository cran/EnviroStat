\name{comb.all}
\alias{comb.all}

\title{
  Enumerate all combinations
}
\description{
  Function to enumerate all combinations of size 'k' from 
  the set {1,.., n}.
}
\usage{
comb.all(n, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    number of elements in the set
  }
  \item{k}{
    number of elements to choose for each combination
  }
}

\value{
  An (n-choose-k) by k matrix with one combination per row.
}

\examples{
comb.all(5, 3)
}
