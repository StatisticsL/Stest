#' @title Zero-thresholding function of quantile affine LASSO
#'
#' @description Get the smallest \eqn{\lambda} value for the quantile affine LASSO estimator to set its penalty to zero.
#'
#' @param y The response vector.
#' @param X The covariate matrix.
#' @param A The linear constraint matrix. If \code{A="diff"} then the total variation test is performed.
#' @param b The constraint vector.
#' @param tau The quantile. The default value is 0.5.
#' @param Rank Logical value. If TRUE, the rankS test is performed; otherwise, the S test is performed. The default value is FALSE.
#' @param outtype If set to 1, returns lambda0. If set to 2, returns lambda0 and indexmax.
#' @param precalculated Either set to TRUE (default value), in which case the precalculatedmatrices function is called,
#' otherwise, precalculated should be the output object of the precalculatedmatrices function.
#' @param q The q/(q-1)-norm of the score used by the test statistic. Either 1 (default value) or 2.
#'
#'
#' @export
#' @return
#' \item{lambda0}{The value of zero-thresholding function.}
#' \item{indexmax}{The index of max value.}
#'
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
#' ztf_aql(y,X,A,b)

ztfaql=function(y,X,A,b,tau=0.5,Rank=F,outtype="1",precalculated=T,q=1){
  n=nrow(X)
  p=ncol(X)

  if(all(A=="diff")){
    temp=reduceXyfromAb(y=y,X=X,A=A,b=b,precalculated=precalculated)
    Xnew=temp$Xnew
    ynew=temp$ynew
    Xout_index=temp$Xout_index
    AAtinv=temp$AAtinv

    obj=rq(ynew~Xnew-1,tau =tau)
    w=obj$dual
    w=-w-tau+1
    beta2hatobj=as.vector(obj$coefficients)
    beta1hatobj=temp$A1inv%*%(b-matrix(c(rep(0,n-2),beta2hatobj),nrow = n-1))
    betahat=rep(NA,p)
    betahat[Xout_index]=beta1hatobj
    betahat[-Xout_index]=beta2hatobj

    res=y-X%*%betahat
    resweighted=(res>0)*tau+(res<0)*(1-tau)
    rankres=rank(abs(res*resweighted))

    if(Rank==T){
      tempv=abs(AAtinv%*%diff(t(X)%*%(w*rankres)))
    }
    else{
      tempv=abs(AAtinv%*%diff(t(X)%*%w))
    }
    indexmax=which.max(tempv)
    if(q==1){
      lambda0=max(tempv)
    }else{
      lambda0=norm(tempv,type = "2")
    }


  }
  else{
    if(!is.matrix(A)){print("A is not a matrix!");break}
    m=nrow(A)
    if(p!=ncol(A)){print("Number of columns of A is different from number of columns of X!");break}
    if(qr(A)$rank != m){print("A is not full row rank!");break}
    if((m>1)|(p>1)){
      temp=reduceXyfromAb(y=y,X=X,A=A,b=b,precalculated=precalculated)
      Xnew=temp$Xnew
      ynew=temp$ynew
      Xout_index=temp$Xout_index
      AAtinvA=temp$AAtinvA

      #############################
      obj=rq(ynew~Xnew-1,tau =tau)
      w=obj$dual
      w=-w-tau+1
      beta2hatobj=as.vector(obj$coefficients)
      beta1hatobj=temp$A1inv%*%(b-temp$A2%*%as.matrix(beta2hatobj,ncol=1))
      betahat=rep(NA,p)
      betahat[Xout_index]=beta1hatobj
      betahat[-Xout_index]=beta2hatobj

      res=y-X%*%betahat
      resweighted=(res>0)*tau+(res<0)*(1-tau)
      rankres=rank(abs(res*resweighted))
      }else {
        w=sign(y-X%*%(b/A))  ###m-p=0
        res=y-X%*%(b/A)
        resweighted=(res>0)*tau+(res<0)*(1-tau)
        rankres=rank(abs(res*resweighted))
      }
    if(Rank==T){
      tempv=abs(AAtinvA%*%t(X)%*%(w*rankres))
    }else{
      tempv=abs(AAtinvA%*%t(X)%*%w)
    }
    indexmax=which.max(tempv)
    if(q==1){
      lambda0=max(tempv)
    }else{
      lambda0=norm(tempv,type = "2")
    }

  }

  out=NULL

  if(outtype=="1"){
    out=lambda0
  }else{
    out$lambda0=lambda0
    out$indexmax=indexmax
  }
  return(out)
}

