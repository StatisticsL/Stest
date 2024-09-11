#' @title Plot for Stest
#'
#' @description Plot the density of lambdas simulated from Monte Carlon, quantile affine LASSO path. And if \code{A="diff"}, also the data, estimated jumps and estimated quantiles.
#'
#' @param objStest The object of Stest.
#' @param main The main for the data. The default value is "Data".
#' @param xlab The xlab for the data and estimated quantiles. The default value is "i".
#' @param ylab The ylab for the data. The default value is "\eqn{y_i}".

#' @export

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
#' obj=Stest(y,X,A,b)
#' plot.Stest(obj)


plot.Stest=function(objStest,main=expression("Data"),xlab=expression("i"),ylab=expression(y[i]),...){
 if(all(objStest$A=="diff")){
   par(mfrow=c(2,2),mai=c(1,1.1,0.4,0.2))
   hist(objStest$outMCH0$lambdas,freq = F,xlim = c(min(objStest$outMCH0$lambdas),max(objStest$outMCH0$lambdas)),
        main=expression(paste("Estimated density of ",S^tau)),xlab = expression(lambda),...)
   axis(1,at=c(objStest$outMCH0$s_critval,objStest$lambda0data),labels = c(expression(c[alpha]^tau),expression(lambda[0]^tau)))

   indexbumps=which(abs(diff(objStest$betahat))>10e-05)
   lindexbumps=length(indexbumps)
   plot(objStest$y,pch="o", main=main,xlab = xlab,ylab = ylab,...)
   for (i in 1:lindexbumps) {
     lines(rep(indexbumps[i]+0.5,2),range(objStest$y),lty=2,lwd=0.5)
   }


   matplot(objStest$lambdas_seq, t(diff(t(objStest$allbetahat))-objStest$b),type="l",lty=1,col="black",xlim = c(min(objStest$outMCH0$lambdas),max(objStest$outMCH0$lambdas)),
           ylab = expression(paste("A",hat(beta),"-b")),xlab = expression(lambda),
           main=expression("Quantile affine LASSO path"),...)#,xaxt="n"
   axis(1,at=c(objStest$outMCH0$s_critval,objStest$lambda0data),labels = c(expression(c[alpha]^tau),expression(lambda[0]^tau)))
   lines(rep(objStest$lambda0data,2),range(t(diff(t(objStest$allbetahat))-objStest$b)), lty=2,col="blue",lwd=2)
   lines(rep(objStest$outMCH0$s_critval,2), range(t(diff(t(objStest$allbetahat))-objStest$b)),lty=2,col="blue",lwd=2)

   plot(objStest$X%*%objStest$betahat,lwd=1,type = "l",ylim=range(objStest$y),xlab =xlab,main = expression("Estimated quantiles"),ylab = expression(paste("X",hat(beta))),...)
 }else{
   par(mfrow=c(2,1),mai=c(1,1.1,0.4,0.2))
   hist(objStest$outMCH0$lambdas,freq = F,xlim = c(min(objStest$outMCH0$lambdas),max(objStest$outMCH0$lambdas)),
        main=expression(paste("Estimated density of ",S^tau)),xlab = expression(lambda),...)
   axis(1,at=c(objStest$outMCH0$s_critval,objStest$lambda0data),labels = c(expression(c[alpha]^tau),expression(lambda[0]^tau)))


   matplot(objStest$lambdas_seq, t(diff(t(objStest$allbetahat))-objStest$b),type="l",lty=1,col="black",xlim = c(min(objStest$outMCH0$lambdas),max(objStest$outMCH0$lambdas)),
           ylab = expression(paste("A",hat(beta),"-b")),xlab = expression(lambda),
           main=expression("quantile affine LASSO path"),...)#,xaxt="n"
   axis(1,at=c(objStest$outMCH0$s_critval,objStest$lambda0data),labels = c(expression(c[alpha]^tau),expression(lambda[0]^tau)))
   lines(rep(objStest$lambda0data,2),range(t(diff(t(objStest$allbetahat))-objStest$b)), lty=2,col="blue",lwd=2)
   lines(rep(objStest$outMCH0$s_critval,2), range(t(diff(t(objStest$allbetahat))-objStest$b)),lty=2,col="blue",lwd=2)
 }
}
