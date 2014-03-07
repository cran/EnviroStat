\name{sinterp}
\alias{sinterp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fitting a thin-plate spline
}
\description{
  Function to compute coefficients for arbitrary dimension thinplate
  spline smooths, that is, smooth mappings from Rq to Rp.
}
\usage{
sinterp(x, y, m = 2, lam = 0, lsq = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An \eqn{nq \times nk}{nq x nk} matrix of knots, where nq = dimension
    of the domain space, and nk = the number of knots. Each column
    represents one knot. That is, the row dimension should be the
    dimension of the domain space, and the column dimension should be
    the number of knot points.
  }
  \item{y}{
    An \eqn{nq \times nk}{nq x nk} matrix of function values at each
    knot point, where np = dimension of the image space. Each data set
    is a column. (Note that the row dimension of y is the number of
    knots which corresponds to the column dimension of x). (As returned
    by \code{\link{Falternate3}}.)
  }
  \item{m}{
    An integer such that \eqn{2*\code{m}} is the order of the spline.
    The default value is 2.
  }
  \item{lam}{
    A vector of real smoothing parameters. If lam is missing or 0,
    sinterp performs interpolation.
  }
  \item{lsq}{
    if \code{TRUE}, subtract least squares fit from \code{y}, returning
    the polynomial coefficients as element \code{b} of the result
  }
}
\value{
  Thin-plate spline solution; more details in Sampson-Guttorp (1992)
  \item{x}{A copy of the \code{x} matrix argument passed to the function.}
  \item{y}{A copy of the \code{y} matrix argument.}
  \item{m}{A copy of the \code{m} argument.}
  \item{lam}{A copy of the \code{lam} argument.}
  \item{lsq}{A copy of the \code{lsq} argument.}
  \item{b}{
    If requested, this contains the coefficients of the least squares
    portion. Each column represents one \code{y} dataset, with the
    coefficients in the same order as in \code{sol}.
  }
  \item{sol}{
    A solution array. The first \eqn{n} components are the coefficients of
    \eqn{u(x-x(i))}. The next \eqn{d+m-1} choose \eqn{m-1} components are the coefficients
    of the interpolating polynomial. These coefficients are ordered by
    increasing order of the total degree of the monomial. Within a group
    of monomials whose total degree is the same, the coefficient of
    monomial \eqn{i} is before the coefficient of monomial \eqn{j} if and only if
    \eqn{x(k)} appears to a higher power in monomial \eqn{i} than it does in
    monomial \eqn{j} for some \eqn{k} such that \eqn{x(k-l)} is not in either monomial for
    any \eqn{l}. Example: if \eqn{m=3} and \eqn{d=4}, the order of the polynomial
    coefficients will be:
    \eqn{1,x1,x2,x3,x4,x1^2,x1*x2,x1*x3,x1*x4,x2^2,x2*x3,x2*x4,x3^2,x3*x4,x4^2}
    The above mess describes each column. The array is \eqn{len \times ny \times nlam}{len x ny x nlam},
    where nlam is the number of lambda values, ny is the number of y
    data vectors, and len is the length of the vector described above.
    In the case that nlam equals 1 (the default), "sol" is returned as
    2-dimensional \eqn{len \times ny}{len x ny} array without the third dimension.
  }
  \item{ainf}{
    A vector of informational integers for the factored "a" matrix.
    Zero's indicate everything is O.K.. If the i-th entry is \eqn{k},
    then the \eqn{k}-th pivot block of the a matrix for the i-th lambda
    value is singular. In the last case, no solution will be computed.
  }
  \item{linf}{
    A vector of informational integers corresponding to the least
    squares solution(if requested). Zero's indicate no problems. If the
    i-th entry is \eqn{k}, then the \eqn{k}-th diagonal entry of the R
    matrix of the QR decomposition of the polynomial matrix is zero.
  }
  \item{f}{
    A copy of the factorial values calculated for various internal
    length determinations.
  }
  \item{a}{
    A copy of the internally generated factored "a" matrix.
  }
}
\references{
  Described in Sampson P.D., and Guttorp, P., "Nonparametric estimation
  of nonstationary spatial covariance structure", in Journal of the
  American Statistical Association, vol 87, pp 108-119, 1992.
}
