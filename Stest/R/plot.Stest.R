#' @title Plot for Stest
#' @description Plot the histogram density of lambdas simulated by Monte Carlo and the quantile affine LASSO path. And if \code{A="diff"}, additionally plots the time series data with estimated jump locations, as well as the estimated regression quantiles.
#' @param objStest The object of Stest.
#' @param main The main for the data. The default value is "Data".
#' @param xlab The xlab for the data and estimated quantiles. The default value is "i".
#' @param ylab The ylab for the data. The default value is "\eqn{y_i}".
#' @param \dots Other graphical parameters to plot
#' @method plot Stest
#' @export
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
#' out1=Stest(y,X,A,b,nlambda=0)
#' out2=Stest(y,X,A,b,nlambda=20)
#' plot(out1)
#' plot(out2)

plot.Stest=function(objStest,main=expression("Data"),xlab=expression("i"),ylab=expression(y[i]),...){
 if(all(objStest$A=="diff")){
   if(length(objStest$lambdas_seq)==0){
     graphics::par(mfrow=c(3,1))
     #####first plot------------------------------------------------
     histresult=hist(objStest$outMCH0$lambdas,plot = F)
     xaxis=c(histresult$breaks,objStest$lambda0data)
     hist(objStest$outMCH0$lambdas,freq = F,xlim = c(min(xaxis),max(xaxis)),xaxt="n",
                    main=expression(paste("Estimated density of ",S^tau)),xlab = expression(lambda),right=FALSE,...)
     axis(1,at=c(seq(min(xaxis),max(xaxis),length=5),objStest$outMCH0$s_critval,objStest$lambda0data),labels = c(seq(min(xaxis),max(xaxis),length=5),expression(c[alpha]^tau),expression(lambda[0]^tau)))

     #####second plot------------------------------------------------
     indexbumps=which(abs(diff(objStest$betahat)-objStest$b)>10e-05)
     lindexbumps=length(indexbumps)
     plot(objStest$y,pch="o", main=main,xlab = xlab,ylab = ylab,...)
     if(lindexbumps>0){
       for (i in 1:lindexbumps) {
         lines(rep(indexbumps[i]+0.5,2),range(objStest$y),lty=2,lwd=0.5)
       }
     }


     #####third plot------------------------------------------------
     plot(objStest$X%*%objStest$betahat,lwd=1,type = "l",ylim=range(objStest$y),xlab =xlab,main = bquote("Estimated quantiles for " * tau ~ "=" ~ .(objStest$tau)),ylab = expression(paste("X",hat(beta))),...)

   }else{
     graphics::par(mfrow=c(2,2),mai=c(1,1.1,0.4,0.2))

     #####first plot------------------------------------------------
     histresult=hist(objStest$outMCH0$lambdas,plot = F)
     xaxis=c(histresult$breaks,objStest$lambda0data)
     hist(objStest$outMCH0$lambdas,freq = F,xlim = c(min(xaxis),max(xaxis)),xaxt="n",xaxt="n",
                    main=expression(paste("Estimated density of ",S^tau)),xlab = expression(lambda),right=FALSE,...)
     axis(1,at=c(seq(min(xaxis),max(xaxis),length=5),objStest$outMCH0$s_critval,objStest$lambda0data),labels = c(seq(min(xaxis),max(xaxis),length=5),expression(c[alpha]^tau),expression(lambda[0]^tau)))

     #####second plot------------------------------------------------
     indexbumps=which(abs(diff(objStest$betahat)-objStest$b)>10e-05)
     lindexbumps=length(indexbumps)
     plot(objStest$y,pch="o", main=main,xlab = xlab,ylab = ylab,...)
     if(lindexbumps>0){
       for (i in 1:lindexbumps) {
         lines(rep(indexbumps[i]+0.5,2),range(objStest$y),lty=2,lwd=0.5)
       }
     }


     #####third plot------------------------------------------------
     homopowerrescalingtemp=objStest$homopowerrescaling
     if(any(is.na(homopowerrescalingtemp))){homopowerrescalingtemp=rep(1,length(objStest$b))}
     matplot(objStest$lambdas_seq, t(diag(homopowerrescalingtemp)%*%(diff(t(objStest$allbetahat)))-objStest$b),type="l",lty=1,col="black",
             ylab = expression(paste("A",hat(beta),"-b")),xlab = expression(lambda),xlim = c(min(xaxis),max(xaxis)),
             main=expression("Quantile affine LASSO path"),xaxt="n",...)#,xaxt="n"
     axis(1,at=c(seq(min(xaxis),max(xaxis),length=5),objStest$outMCH0$s_critval,objStest$lambda0data),labels = c(seq(min(xaxis),max(xaxis),length=5),expression(c[alpha]^tau),expression(lambda[0]^tau)))
     lines(rep(objStest$lambda0data,2),range(t(diag(homopowerrescalingtemp)%*%(diff(t(objStest$allbetahat)))-objStest$b)), lty=2,col="blue",lwd=2)
     lines(rep(objStest$outMCH0$s_critval,2), range(t(diag(homopowerrescalingtemp)%*%(diff(t(objStest$allbetahat)))-objStest$b)),lty=2,col="blue",lwd=2)

     #####fourth plot------------------------------------------------
     plot(objStest$X%*%objStest$betahat,lwd=1,type = "l",ylim=range(objStest$y),xlab =xlab,main =
            bquote("Estimated quantiles for " * tau ~ "=" ~ .(objStest$tau)),ylab = expression(paste("X",hat(beta))),...)
   }
 }else{
   if(length(objStest$lambdas_seq)==0){
     graphics::par(mfrow=c(1,1))

     #####first plot------------------------------------------------
     histresult=hist(objStest$outMCH0$lambdas,plot = F)
     xaxis=c(histresult$breaks,objStest$lambda0data)
     hist(objStest$outMCH0$lambdas,freq = F,xlim = c(min(xaxis),max(xaxis)),xaxt="n",
          main=expression(paste("Estimated density of ",S^tau)),xlab = expression(lambda),right=FALSE,...)
     axis(1,at=c(seq(min(xaxis),max(xaxis),length=5),objStest$outMCH0$s_critval,objStest$lambda0data),labels = c(seq(min(xaxis),max(xaxis),length=5),expression(c[alpha]^tau),expression(lambda[0]^tau)))

   }else{
     graphics::par(mfrow=c(2,1),mai=c(1,1.1,0.4,0.2))

     #####first plot------------------------------------------------
     histresult=hist(objStest$outMCH0$lambdas,plot = F)
     xaxis=c(histresult$breaks,objStest$lambda0data)
     hist(objStest$outMCH0$lambdas,freq = F,xlim = c(min(xaxis),max(xaxis)),xaxt="n",
          main=expression(paste("Estimated density of ",S^tau)),xlab = expression(lambda),right=FALSE,...)
     axis(1,at=c(seq(min(xaxis),max(xaxis),length=5),objStest$outMCH0$s_critval,objStest$lambda0data),labels = c(seq(min(xaxis),max(xaxis),length=5),expression(c[alpha]^tau),expression(lambda[0]^tau)))


     #####second plot------------------------------------------------
     homopowerrescalingtemp=objStest$homopowerrescaling
     if(any(is.na(homopowerrescalingtemp))){homopowerrescalingtemp=rep(1,length(objStest$b))}
     colabeta=dim(objStest$A%*%t(objStest$allbetahat))[2]
     matplot(objStest$lambdas_seq, t(diag(homopowerrescalingtemp)%*%objStest$A%*%t(objStest$allbetahat)-matrix(rep(objStest$b,colabeta),ncol = colabeta)) ,type="l",lty=1,col="black",xlim = c(min(xaxis),max(xaxis)),
             ylab = expression(paste("A",hat(beta),"-b")),xlab = expression(lambda),
             main=expression("Quantile affine LASSO path"),xaxt="n",...)
     axis(1,at=c(seq(min(xaxis),max(xaxis),length=5),objStest$outMCH0$s_critval,objStest$lambda0data),labels = c(seq(min(xaxis),max(xaxis),length=5),expression(c[alpha]^tau),expression(lambda[0]^tau)))
     lines(rep(objStest$lambda0data,2),range(t(diag(homopowerrescalingtemp)%*%objStest$A%*%t(objStest$allbetahat)-matrix(rep(objStest$b,colabeta),ncol = colabeta))), lty=2,col="blue",lwd=2)
     lines(rep(objStest$outMCH0$s_critval,2), range(t(diag(homopowerrescalingtemp)%*%objStest$A%*%t(objStest$allbetahat)-matrix(rep(objStest$b,colabeta),ncol = colabeta))),lty=2,col="blue",lwd=2)
   }
 }
}
