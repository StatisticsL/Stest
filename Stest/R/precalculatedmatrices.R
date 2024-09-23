#' @title Precalculated matrices to reduce the constrained regression
#'
#' @description For the linear constraint \eqn{A\beta=b},
#' the matrix A must be decomposed into A=[A1, A2], where A1 is nonsingular.
#' Related matrices can also be precalculated to avoid repeating this operation several times.
#' This option is particularly useful for experienced users.
#
#'
#'
#' @param A The linear constraint matrix. If \code{A="diff"} then the total variation test is performed.
#' @param n The sample size. The default value is NA. If \code{A="diff"}, then the sample size n must be provided.

#' @export
#' @return
#' \item{A1}{The nonsingular matrix A1}
#' \item{A2}{The A2 matrix such that A=[A1, A2].}
#' \item{Xout_index}{The indexes of columns of A leading to a nonsingular matrix A1.}
#' \item{A1inv}{The inverse of A1.}
#' \item{AAtinv}{The inverse of \eqn{AA^{T}}.}
#' \item{AAtinvA}{\eqn{(AA^{T})^{-1}A}}
#'
#' @examples
#' p=20
#' A=diag(5)
#' A=cbind(A,matrix(0,5,p-5))
#' A[1,2]=A[2,3]=A[3,4]=A[4,5]=A[5,6]=-1
#' out=precalculatedmatrices(A)
#' ##############
#' A="diff"
#' n=10
#' out=precalculatedmatrices(A,n)


precalculatedmatrices<-function(A,n=NA){
  if(all(A=="diff")){
    if(is.na(n)){stop("You should provide the sample size!")}
    Xout_index=seq(1,n-1,by=1)
    A=cbind(-1*diag(n-1),rep(0,n-1))+cbind(rep(0,n-1),diag(n-1))
    A1=A[,Xout_index]
    A2=A[,-Xout_index]
    A1inv=upper.tri(diag(n-1),diag = T)
    A1inv[which(A1inv==T)]=-1
    AAtinv=matrix(NA,n-1,n-1)
    for (i in 1:(n-1)) {
      for (j in 1:(n-1)) {
        AAtinv[i,j]=min(i,j) -i*j/n
      }
    }
    AAtinvA=NA
  }else{
    out=gaussianElimination(A)
    Xout_index=(which(apply(out,2,sum)==1))
    A1=A[,Xout_index]
    A2=A[,-Xout_index]
    A1inv=solve(A1)
    AAtinv=solve(A%*%t(A))
    AAtinvA=AAtinv%*%A
  }
  out=NULL
  out$Xout_index=Xout_index
  out$A1=A1
  out$A2=A2
  out$A1inv=A1inv
  out$AAtinv=AAtinv
  out$AAtinvA=AAtinvA

  return(out)
}
