#' @title Monte Carlo estimation of the critial value of the S test.
#'
#' @description Estimate the \eqn{c_{\alpha}^{\tau}} value by Monte Carlo by simulated the distribution of the test statistic S under H0.
#'
#' @param X The covariate matrix.
#' @param A The linear constraint matrix. If \code{A="diff"} then the total variation test is performed.
#' @param b The constraint vector.
#' @param tau The quantile. The default value is 0.5.
#' @param alpha The test level. The default value is 0.05.
#' @param M The number of Monte Carlo runs. The default value is 1000.
#' @param Rank Logical value. If TRUE, the rankS test is performed; otherwise, the S test is performed. The default value is FALSE.
#' @param precalculated Either set to TRUE (default value), in which case the precalculatedmatrices function is called,
#' otherwise, precalculated should be the output object of the precalculatedmatrices function.
#' @param q The q/(q-1)-norm of the score used by the test statistic. Either 1 (default value corresponding to a LASSO penalty) or 2 (corresponding to a group LASSO penalty).
#'
#' @export
#' @return
#' \item{lambdas}{All M samples of the test statistic S under H0 simulated by Monte Carlo.}
#' \item{s_critval}{The estimated critical value \eqn{c_{\alpha}^{\tau}} of the test by taking the upper \eqn{(1-\alpha)}-quantile of the M samples of the test statistic S under H0.}
#'
#' @examples
#' n=100
#' p=20
#' X=matrix(rt(n*p, df=2),n,p)
#' beta=stats::rnorm(p)
#' A=diag(5)
#' A=cbind(A,matrix(0,5,p-5))
#' A[1,2]=A[2,3]=A[3,4]=A[4,5]=A[5,6]=-1
#' b=A%*%beta
#' MCH0(X,A,b)

MCH0<-function(X,A,b,tau=0.5,alpha=0.05,M=1000,Rank=F,precalculated=T,q=1){
  n=nrow(X)
  p=ncol(X)
  if(isTRUE(precalculated)){
    precalculated=precalculatedmatrices(A=A,n=n)
  }
  Xout_index=precalculated$Xout_index
  A1=precalculated$A1
  A2=precalculated$A2
  A1inv=precalculated$A1inv
  AAtinv=precalculated$AAtinv
  AAtinvA=precalculated$AAtinvA

  beta=rep(0,p)
  beta[Xout_index]=A1inv%*%b

  Y=matrix(stats::rnorm(n*M),n,M)+as.vector(X%*%beta)
  lambdas=rep(NA,M)
  for (j in 1:M) {lambdas[j]=ztfaql(y=Y[,j],X=X,A=A,b=b,tau=tau,Rank=Rank,outtype="1",precalculated=precalculated,q=q)}
  s_critval=stats::quantile(lambdas,1-alpha)
  out=NULL
  out$lambdas=lambdas
  out$s_critval=s_critval
  return(out)
}
