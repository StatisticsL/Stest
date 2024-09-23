#' @title Rescale the \code{W}
#'
#' @description Rescale \code{W} using homopowerrescaling.
#'
#' @param X The covariate matrix.
#' @param A The linear constraint matrix. If \code{A="diff"} then the total variation test is performed.
#' @param b The constraint vector.
#' @param homopowerrescaling Weights used to rescale the rows of A and b.
#' @param tau The quantile. The default value is 0.5.
#' @param alpha The quantile \eqn{(1-\alpha)} to rescale. The default value is 0.05.
#' @param M The number of Monte Carlo runs. The default value is 1000.
#' @param Rank Logical value. If TRUE, the rankS test is p erformed; otherwise, the S test is performed. The default value is FALSE.
#' @param precalculated Either set to TRUE (default value), in which case the precalculatedmatrices function is called,
#' otherwise, precalculated should be the output object of the precalculatedmatrices function.
#'
#' @export
#'
#' @return
#' \item{allwvec}{Matrix \code{W} simulated by Monte Carlo used for rescaling A and b.}
#' \item{allwvecrescaled}{The rescaled vectors \code{W} before calculated the sampled test statistic under \eqn{H_0}.}
#' \item{homopowerrescaling}{Weights used to rescale the rows of A and b.}
#' @examples
#' n=100
#' p=20
#' X=matrix(rnorm(n*p),n,p)
#' beta=rnorm(p)
#' A=diag(5)
#' A=cbind(A,matrix(0,5,p-5))
#' A[1,2]=A[2,3]=A[3,4]=A[4,5]=A[5,6]=-1
#' b=A%*%beta
#' rescaleW(X,A,b)

rescaleW<-function(X,A,b,homopowerrescaling=NA,tau=0.5,alpha=0.05,M=1000,Rank=F,precalculated=T){
  n=nrow(X)
  p=ncol(X)
  if(isTRUE(precalculated)){
    precalculated=precalculatedmatrices(A=A,n=n)
  }
  if(any(is.na(homopowerrescaling))){
    homopowerrescaling=rescaleAb(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=Rank,precalculated=precalculated)$homopowerrescaling
  }
  Xout_index=precalculated$Xout_index
  A1=precalculated$A1
  A2=precalculated$A2
  A1inv=precalculated$A1inv
  AAtinv=precalculated$AAtinv
  AAtinvA=precalculated$AAtinvA

  beta=rep(0,p)
  beta[Xout_index]=A1inv%*%b

  if(n!=p){
    Y=matrix(rnorm(n*M),n,M)+as.vector(X%*%beta)
  }else{
    if(all(X==diag(n))){
      Y=matrix(rnorm(n*M),n,M)+as.vector(beta)
    }else{
      Y=matrix(rnorm(n*M),n,M)+as.vector(X%*%beta)
    }
  }

  allwvec=matrix(NA, M, nrow(A1))
  for (j in 1:M) {
    out=ztfaql(y=Y[,j],X=X,A=A,b=b,tau=tau,Rank=Rank,outtype="2",precalculated=precalculated,q=1)
    allwvec[j,]=out$wvec
  }

  allwvecrescaled=t(diag(1/homopowerrescaling)%*%t(allwvec))

  out=NULL
  out$allwvec=allwvec
  out$allwvecrescaled=allwvecrescaled
  out$homopowerrescaling=homopowerrescaling
  return(out)
}
