# install_github("StatisticsL/Stest/Stest")
library(Stest)
###----------------------------------------------
####F_test function-----------------------
F_test<-function(y,X,A,b,alpha=0.05,AXtXinvAtinv=NA){
  n=nrow(X)
  p=ncol(X)
  m=nrow(A)
  betahat=as.vector(lm(y~X-1)$coefficients)

  if(!is.matrix(AXtXinvAtinv)){AXtXinvAtinv=solve(A%*%solve(t(X)%*%X)%*%t(A))}
  diff=A%*%betahat-b
  s2=sum((y-X%*%betahat)^2)/(n-p)
  Fdata=(t(diff)%*%AXtXinvAtinv%*%diff)/(m*s2)
  F_critval=qf(1-alpha,m,n-p)

  out=NULL
  out$p_value=1-pf(Fdata,m,n-p)
  out$test_value=as.numeric(Fdata>F_critval)

  return(out)
}


###----------------------------------------------
n=100
p=20
X=matrix(rt(n*p, df=2),n,p)
beta=rnorm(p)
A=diag(5)
A=cbind(A,matrix(0,5,p-5))
A[1,2]=A[2,3]=A[3,4]=A[4,5]=A[5,6]=-1
b=A%*%beta
alpha=0.01
tau=0.5
M=10000
#######################################################
precalculated=precalculatedmatrices(A=A,n=n)
out_rescaled=rescaleAb(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,precalculated=precalculated)
homopowerrescaling=out_rescaled$homopowerrescaling
rankout_rescaled=rescaleAb(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=T,precalculated=precalculated)
rankhomopowerrescaling=rankout_rescaled$homopowerrescaling
Sout_MCH0=MCH0(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,rescale=T,precalculated = precalculated,homopowerrescaling = homopowerrescaling,q=1)
rankSout_MCH0=MCH0(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=T,rescale=T,precalculated = precalculated,homopowerrescaling = rankhomopowerrescaling,q=1)

AXtXinvAtinv=solve(A%*%solve(t(X)%*%X)%*%t(A))

numiter=10000
dfs=c(1,2,3,4,5,10,100); ldfs=length(dfs)
resultsS=matrix(NA,ldfs,numiter)
resultsrankS=matrix(NA,ldfs,numiter)
resultsF=matrix(NA,ldfs,numiter)

for(j in 1:ldfs){
  print(j)
  for (iter in 1:numiter) {
    err=rt(n,dfs[j])
    y=X%*%beta+err

    resultsS[j,iter]=Stest(y=y,X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,rescale=T,outMCH0=Sout_MCH0, precalculated = precalculated,homopowerrescaling = homopowerrescaling,q=1,nlambda = 0)$test_value
    resultsF[j,iter]=F_test(y=y,X=X,A=A,b=b,alpha=alpha,AXtXinvAtinv=AXtXinvAtinv)$test_value
    resultsrankS[j,iter]=Stest(y=y,X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=T,rescale = T,outMCH0=rankSout_MCH0,precalculated = precalculated,homopowerrescaling = rankhomopowerrescaling,q=1,nlambda = 0)$test_value

  }
}


####plot----------------------------------------------------
plot(apply(resultsS,1,mean),ylab="Level", xlab=expression("df"),ylim=c(0,0.1),xaxt="n",main = expression("Levels of three tests"),type="b",pch=19)
axis(1,at=seq(1,ldfs,by=1),labels = dfs)
lines(apply(resultsrankS,1,mean),pch=3,col="red")
lines(apply(resultsF,1,mean), pch=5,col="blue")
abline(c(alpha, 0), lty=2)
legend("topleft",legend = c(expression(paste(infinity,"-S test")), expression(paste(infinity,"-rankS test")),"F-test"),col=c("black","red","blue"),lty=c(1,1,1),pch = c(19,NA,NA))

save(precalculated,out_rescaled,homopowerrescaling,rankout_rescaled,rankhomopowerrescaling,Sout_MCH0,rankSout_MCH0,
    X,beta,resultsS,resultsrankS,resultsF,file="levelrobustMC.RData")
