#' @title Calcultion of regression quantile affine LASSO estimate.
#'
#' @description For a given penalty parameter lambda, calculate the solution of the regression quantile affine LASSO optimization prob.
#'lem.
#' @param y The response vector.
#' @param X The covariate matrix.
#' @param A The linear constraint matrix. If \code{A="diff"} then the total variation test is performed.
#' @param b The constraint vector.
#' @param tau The quantile. The default value is 0.5.
#' @param lambda The penalty parameter.
#' @param q The q/(q-1)-norm of the score used by the test statistic. Either 1 (default value corresponding to a LASSO penalty) or 2 (corresponding to a group LASSO penalty).
#'
#' @export
#' @return
#' \item{loss_value}{The loss value of regression quantile affine LASSO.}
#' \item{betahat}{The solution of regression quantile affine LASSO.}
#'
#' @examples
#' n=100
#' p=20
#' X=matrix(rt(n*p, df=2),n,p)/sqrt(n)
#' beta=c(c(3,3,4,4),rnorm(p-4,sd=5))
#' A=diag(5)
#' A=cbind(A,matrix(0,5,p-5))
#' A[1,2]=A[2,3]=A[3,4]=A[4,5]=A[5,6]=-1
#' b=A%*%beta
#' y=X%*%beta+rt(n,1)
#' rqaffineLASSO(y,X,A,b,lambda=1)


rqaffineLASSO <- function(y,X,A,b,tau=0.5,lambda,q=1){
  p=ncol(X)
  n=nrow(X)
  if(all(A=="diff")){A=cbind(-1*diag(n-1),rep(0,n-1))+cbind(rep(0,n-1),diag(n-1))}
  betaHat <- CVXR::Variable(p)
  quant_loss <- function(u, tau) { 0.5 * abs(u) + (tau - 0.5) * u }
  ###
  if(q==1){
    objective <- sum(quant_loss(y-X%*%betaHat, tau)) + lambda*sum(abs(A%*%betaHat-b) )
  }else{
    objective <- sum(quant_loss(y-X%*%betaHat, tau)) + lambda*norm((A%*%betaHat-b),type = "2")
  }
  problem <- CVXR::Problem(CVXR::Minimize(objective))
  result <- CVXR::solve(problem, adaptive_rho = FALSE, verbose = F,num_iter=10e+6)
  out=NULL
  out$loss_value=result$value
  out$betahat=result$getValue(betaHat)
  return(out)
}
