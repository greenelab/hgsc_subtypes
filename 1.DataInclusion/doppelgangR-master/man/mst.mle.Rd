\name{mst.mle}
\alias{mst.mle}
\alias{st.mle}
\title{
Maximum likelihood estimation for a (multivariate) skew-t distribution
}
\description{
Fits a skew-t (ST) or multivariate skew-t (MST) distribution to data, or
fits a linear regression model with (multivariate) skew-t errors,
using maximum likelihood estimation.  Functions copied from \code{sn}
CRAN library v0.4.18 because they were later deprecated in that library.
}
\usage{
mst.mle(X, y, freq, start, fixed.df=NA, trace=FALSE,
   algorithm = c("nlminb","Nelder-Mead", "BFGS", "CG", "SANN"), control=list())
st.mle(X, y, freq, start, fixed.df=NA, trace=FALSE,
   algorithm = c("nlminb","Nelder-Mead", "BFGS", "CG", "SANN"), control=list())

}
\arguments{
\item{y}{
a matrix (for \code{mst.mle}) or a vector (for \code{st.mle}).
If \code{y} is a matrix, rows refer to observations, and columns to
components of the multivariate distribution.
}
\item{X}{
  a matrix of covariate values.
  If missing, a one-column matrix of 1's is created; otherwise,
  it must have the same number of rows of \code{y}.
  If \code{X} is supplied, then it must include a column of 1's.
}
\item{freq}{
a vector of weights.
If missing, a vector of 1's is created; otherwise
it must have  length equal to the number of rows of \code{y}.
}
\item{start}{
  for \code{mst.mle}, a list contaning the components
  \code{beta},\code{Omega},  \code{alpha},
\code{df} of the type described below; for \code{st.mle}, a vector whose
components contain analogous ingredients as before, with the exception
that the scale parameter is the square root of \code{Omega}.  In both
cases, the \code{dp} component of the returned list from a previous call
has the required format and it can be used as a new \code{start}. If the
\code{start} parameter is missing, initial values are selected by the
function.

}
\item{fixed.df}{
a scalar value containing the degrees of freedom (df), if these must
be taked as fixed, or \code{NA} (default value) if \code{df} is a parameter
to be estimated.
}
\item{trace}{
logical value which controls printing of the algorithm convergence.
If \code{trace=TRUE}, details are printed. Default value is \code{FALSE}.
}
\item{algorithm}{
a character string which selects the numerical optimization procedure
used to maximize the loglikelihood function. If this string is set
equal to \code{"nlminb"}, then this function is called; in all other cases,
\code{optim} is called, with \code{method} set equal to the given string.
Default value is \code{"nlminb"}.
}
\item{control}{
  this parameter is passed to the chose optimizer, either \code{nlminb}
  or  \code{optim}; see the documentation of this function for its usage.
}}
\value{
A list containing the following components:

\item{call}{
a string containing the calling statement.
}
\item{dp}{
  for \code{mst.mle}, this is a list containing the direct parameters
  \code{beta}, \code{Omega}, \code{alpha}.
Here, \code{beta} is a matrix of regression coefficients with
\code{dim(beta)=c(ncol(X),ncol(y))}, \code{Omega} is a covariance matrix
of order \code{ncol(y)}, \code{alpha} is a vector of shape parameters of
length \code{ncol(y)}.  For \code{st.mle}, \code{dp} is a vector of
length \code{ncol(X)+3}, containing \code{c(beta, omega, alpha, df)}, where
\code{omega} is the square root of \code{Omega}.
}
\item{se}{
a list containing the components \code{beta}, \code{alpha}, \code{info}.
Here, \code{beta} and \code{alpha} are the standard errors for the
corresponding point estimates;
\code{info} is the observed information matrix for the working parameter,
as explained below.
}
\item{algorithm}{
  the list returned by the chose optimizer, either \code{nlminb}
  or  \code{optim}, plus an item with the \code{name} of the selected
  algorithm; see the documentation of either \code{nlminb}
  or  \code{optim}  for explanation  of the other components.
}}
\details{
If \code{y} is a vector and it is supplied to \code{mst.mle}, then
it is converted to a one-column matrix, and a scalar skew-t distribution
is fitted. This is also the mechanism used by \code{st.mle}
which is simply an interface to \code{mst.mle}.

The parameter \code{freq} is intended for use with grouped data,
setting the values of \code{y} equal to the central values of the
cells; in this case the resulting estimate is an approximation
to the exact maximum likelihood estimate. If \code{freq} is not
set, exact maximum likelihood estimation is performed.

% To fit a scalar skew-t distribution to grouped data by exact
% maximum likelihood estimation, use \code{st.mle.grouped}.

Numerical search of the maximum likelihood estimates is performed in a
suitable re-parameterization of the original parameters with aid of the
selected optimizer (\code{nlminb} or \code{optim}) which is supplied
with the derivatives of the log-likelihood function. Notice that, in
case the optimizer is \code{optim}), the gradient may or may not be
used, depending on which specific method has been selected.  On exit
from the optimizer, an inverse transformation of the parameters is
performed. For a specific description on the re-parametrization adopted,
see Section 5.1 and Appendix B of Azzalini \& Capitanio (2003).

}
\section{Background}{
The family of multivariate skew-t distributions is an extension of the
multivariate Student's t family, via the introduction of a \code{shape}
parameter which regulates skewness; when \code{shape=0}, the skew-t
distribution reduces to the usual t distribution.
When \code{df=Inf} the distribution reduces to the multivariate skew-normal
one; see \code{dmsn}. See the reference below for additional information.
}
\references{
Azzalini, A. and Capitanio, A. (2003).
  Distributions generated by perturbation of symmetry
  with emphasis on a multivariate skew \emph{t} distribution.
  The full version of the paper published in abriged form in
  \emph{J.Roy. Statist. Soc. B} \bold{65}, 367--389,
  is available at \url{http://azzalini.stat.unipd.it/SN/se-ext.ps}
}
\seealso{
  \code{\link{dst}}
}
\examples{
dat <- rt(100, df=5, ncp=100)
fit <- st.mle(y=dat)
fit
}
\keyword{distribution}
\keyword{regression}

