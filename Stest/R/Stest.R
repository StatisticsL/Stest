#' @title S test
#'
#' @description Calculate the test function value (0 or 1) and p-value for the S test applied to the data for the desired linear null hypopthesis.
#'
#' @param y The response vector.
#' @param X The covariate matrix.
#' @param A The linear constraint matrix. If \code{A="diff"} then the total variation test is performed.
#' @param b The constraint vector.
#' @param tau The quantile. The default value is 0.5.
#' @param alpha The test level. The default value is 0.05.
#' @param M The number of Monte Carlo runs. The default value is 1000.
#' @param Rank Logical value. If TRUE, the rankS test is performed; otherwise, the S test is performed. The default value is FALSE.
#' @param outMCH0 The Monte Carlo experiments outcome providing \eqn{c_{\alpha}^{\tau}} and all M sampled test statistics value. The default value is NA.
#' @param precalculated Either set to TRUE (default value), in which case the precalculatedmatrices function is called,
#' otherwise, precalculated should be the output object of the precalculatedmatrices function.
#' @param q The q/(q-1)-norm of the score used by the test statistic. Either 1 (default value) or 2.
#' @param nlambda The number of lambda values to compute the path. The default is 20. Set it to 0 if not need the path.
#' @export
#' @return
#' \item{p_value}{The p value of S test.}
#' \item{test_value}{The test value of S test.}
#' \item{precalculated}{The object of the precalculatedmatrices function.}
#' \item{outMCH0}{\eqn{c_{\alpha}^{\tau}} and all M sampled test statistics value.}
#' \item{lambdas_seq}{Lambdas for the quantile affine LASSO path.}
#' \item{allbetahat}{All betahat corresponding to the lambdas_seq.}
#' \item{betahat}{The betahat corresponding to the critial value.}
#' \item{lambda0data}{The lambda from the data.}
#' \item{X}{The covariate matrix.}
#' \item{y}{The response vector.}
#' \item{A}{The linear constraint matrix.}
#' \item{b}{The constraint vector.}
#' \item{tau}{The quantile.}
#' \item{alpha}{The test level.}
#' \item{q}{The q/(q-1)-norm of the score used by the test statistic.}
#' @examples
#' n=100
#' p=20
#' X=matrix(rt(n*p, df=2),n,p)
#' beta=rnorm(p)
#' A=diag(5)
#' A=cbind(A,matrix(0,5,p-5))
#' A[1,2]=A[2,3]=A[3,4]=A[4,5]=A[5,6]=-1
#' b=A%*%beta
#' y=X%*%beta+rt(n,1)
#' Stest(y,X,A,b,nlambda=0)

Stest<-function(y,X,A,b,tau=0.5,alpha=0.05,M=1000,Rank=F,outMCH0=NA,precalculated=T,q=1,nlambda=20){
  n=length(y)
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

  if(any(is.na(outMCH0))){
    outMCH0=MCH0(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=Rank,precalculated=precalculated,q=q)
  }
  lambdas=outMCH0$lambdas
  s_critval=outMCH0$s_critval

  lambda0data=ztfaql(y=y,X=X,A=A,b=b,tau=tau,Rank=Rank,outtype="1",precalculated=precalculated,q=q)

  p_value=mean(lambda0data<lambdas)
  test_value=as.numeric(lambda0data>s_critval)

  lambdas_seq=seq(min(lambdas),max(lambdas),length=nlambda)
  if(nlambda>0){
    allbetahat=matrix(NA, nrow=nlambda,p)
    for(i in 1:nlambda){
      outi=rqaffineLASSO(y=y,X=X,A=A,b=b,tau=tau,lambdas_seq[i],q=q)
      allbetahat[i,]=outi$betahat
    }
    outi=rqaffineLASSO(y=y,X=X,A=A,b=b,tau=tau,lambda = s_critval,q=q)
    betahat=outi$betahat
  }else{
    allbetahat=NA
    betahat=NA
  }



  #########return---------------------------------
  out=NULL
  out$p_value=p_value
  out$test_value=test_value
  out$precalculated=precalculated
  out$outMCH0=outMCH0
  out$lambdas_seq=lambdas_seq
  out$allbetahat=allbetahat
  out$betahat=betahat
  out$lambda0data=lambda0data
  out$X=X
  out$y=y
  out$A=A
  out$b=b
  out$tau=tau
  out$alpha=alpha
  out$q=q

  return(out)
}
