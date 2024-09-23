#' @title Transform a LAD optimization constrained to \eqn{A\beta=b} to an unconstrained LAD of lower dimension
#'
#' @description Returns the new y and new X via Gaussian elimitation of the equivalent LAD optimization without constraints.
#'
#' @param y The response vector.
#' @param X The covariate matrix.
#' @param A The linear constraint matrix. If \code{A="diff"} then the total variation test is performed.
#' @param b The constraint vector.
#' @param precalculated Either set to TRUE (default value), in which case the precalculatedmatrices function is called,
#' otherwise, precalculated should be the output object of the precalculatedmatrices function.
#'
#' @export
#' @return
#' \item{ynew}{The new y.}
#' \item{Xnew}{The new X.}
#' \item{A1}{The nonsingular matrix A1}
#' \item{A2}{The A2 matrix such that A=[A1, A2].}
#' \item{Xout_index}{The indexes of columns of A leading to a nonsingular matrix A1.}
#' \item{A1inv}{The inverse of A1.}
#' \item{AAtinv}{The inverse of \eqn{AA^{T}}.}
#' \item{AAtinvA}{\eqn{(AA^{T})^{-1}A}}
#'
#' @examples
#' n=100
#' p=20
#' X=matrix(rnorm(n*p),n,p)
#' beta=rnorm(p)
#' A=diag(5)
#' A=cbind(A,matrix(0,5,p-5))
#' A[1,2]=A[2,3]=A[3,4]=A[4,5]=A[5,6]=-1
#' b=A%*%beta
#' y=X%*%beta+rnorm(n,1)
#' reduceXyfromAb(y,X,A,b)

reduceXyfromAb=function(y,X,A,b,precalculated=T){
  n=length(y)
  if(isTRUE(precalculated)){
    precalculated=precalculatedmatrices(A=A,n=n)
  }
  Xout_index=precalculated$Xout_index
  A1=precalculated$A1
  A2=precalculated$A2
  A1inv=precalculated$A1inv
  AAtinv=precalculated$AAtinv
  AAtinvA=precalculated$AAtinvA


  prodX1A1inv=X[,Xout_index]%*%A1inv
  ###########################
  ynew=y-as.vector(prodX1A1inv%*%b)
  Xnew=X[,-Xout_index]-prodX1A1inv%*%A2

  out=NULL
  out$ynew=ynew
  out$Xnew=Xnew
  out$A1inv=A1inv
  out$A2=A2
  out$A1=A1
  out$AAtinv=AAtinv
  out$AAtinvA=AAtinvA
  out$Xout_index=Xout_index
  return(out)
}
