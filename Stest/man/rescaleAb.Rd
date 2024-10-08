% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rescaleAb.R
\name{rescaleAb}
\alias{rescaleAb}
\title{Monte Carlo estimation of the matrix to rescale A and b}
\usage{
rescaleAb(
  X,
  A,
  b,
  tau = 0.5,
  alpha = 0.05,
  M = 1000,
  Rank = F,
  precalculated = T
)
}
\arguments{
\item{X}{The covariate matrix.}

\item{A}{The linear constraint matrix. If \code{A="diff"} then the total variation test is performed.}

\item{b}{The constraint vector.}

\item{tau}{The quantile. The default value is 0.5.}

\item{alpha}{The quantile \eqn{(1-\alpha)} to rescale. The default value is 0.05.}

\item{M}{The number of Monte Carlo runs. The default value is 1000.}

\item{Rank}{Logical value. If TRUE, the rankS test is p erformed; otherwise, the S test is performed. The default value is FALSE.}

\item{precalculated}{Either set to TRUE (default value), in which case the precalculatedmatrices function is called,
otherwise, precalculated should be the output object of the precalculatedmatrices function.}
}
\value{
\item{rescaledA}{The rescaled \code{A}.}
\item{rescaledb}{The rescaled \code{b}.}
\item{homopowerrescaling}{Weights used to rescale the rows of A and b.}
}
\description{
Estimate the homopowerrescaling matrix.
}
\examples{
n=100
p=20
X=matrix(rnorm(n*p),n,p)
beta=rnorm(p)
A=diag(5)
A=cbind(A,matrix(0,5,p-5))
A[1,2]=A[2,3]=A[3,4]=A[4,5]=A[5,6]=-1
b=A\%*\%beta
rescaleAb(X,A,b)

}
