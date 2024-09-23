#' @title Summary for Stest
#'
#' @description Reports the estimated p-value and the decision of the S Stest .
#'
#' @param objStest The object of Stest.

#' @export

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
#' out=Stest(y,X,A,b)
#' summaryStest(out)

summaryStest<-function(objStest){
  print(paste("The estimated p-value by Monte Carlo is ",objStest$p_value,sep = ""))
  if(objStest$test_value==0){
    print(paste("The null hypothesis was not rejected at the level alpha=",objStest$alpha,sep=""))
  }else{
    print(paste("The null hypothesis was rejected at the level alpha=",objStest$alpha,sep=""))
  }
}
