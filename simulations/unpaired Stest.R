# install_github("StatisticsL/Stest/Stest")
library(Stest)
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


###-----------------------------------------------------------------------
na=20
nb=20
p=2
A=matrix(c(0,1),1,2)
b=0
mu=2
X=matrix(c(rep(1,na+nb),rep(0,na),rep(1,nb)),(na+nb),2,byrow = F)
alpha=0.05
tau=0.5
M=10000
#######################################################
precalculated=precalculatedmatrices(A=A)
Sout_MCH0=MCH0(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,rescale=F,precalculated = precalculated,q=1)
rankSout_MCH0=MCH0(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=T,rescale=F,precalculated = precalculated,q=1)

##F-test---------------------
AXtXinvAtinv=solve(A%*%solve(t(X)%*%X)%*%t(A))

##Roger's test-----------------
R=t(A)
V=matrix(0,2,1)
V[1,1]=1
p1=1
p2=1
X1=X[,1:p1]
X2=X[,-(1:p1)]
X1new=X%*%V%*%solve(t(V)%*%V)
X2new=X%*%R%*%solve(t(R)%*%R)
yout=X%*%R%*%solve(t(R)%*%R)%*%b
###-------------------------------


numiter=10000

deltas=seq(0,1,by=0.1)
ldeltas=length(deltas)
resultsS=matrix(NA,ldeltas,numiter)
resultsF=matrix(NA,ldeltas,numiter)
resultsrankS=matrix(NA,ldeltas,numiter)
resultsRoger=matrix(NA,ldeltas,numiter)

for (j in 1:ldeltas) {
  print(j)
  for (iter in 1:numiter) {
    erra=rt(na,df=3)
    errb=rt(nb,df=3)
    aa=mu+erra
    delta=deltas[j]
    bb=mu+delta+errb

    y=c(aa,bb)
    groupg=c(rep(1,na),rep(2,nb))

    resultsS[j,iter]=Stest(y=y,X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,rescale=F,outMCH0=Sout_MCH0, precalculated = precalculated,q=1,nlambda = 0)$test_value
    resultsF[j,iter]=F_test(y=y,X=X,A=A,b=b,alpha=alpha,AXtXinvAtinv=AXtXinvAtinv)$test_value
    resultsrankS[j,iter]=Stest(y=y,X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=T,rescale=F,outMCH0=rankSout_MCH0,precalculated = precalculated,q=1,nlambda = 0)$test_value

    ##Roger's test----
    Ynew=y-yout
    fit0=rq(Ynew~X1new-1)
    fit1=rq(Ynew~X1new+X2new-1)
    out=anova(fit1,fit0,test = "rank")
    resultsRoger[j,iter]=sum((out$table[3]>=qchisq(1-alpha, df=as.numeric(out$table[1]))))
    ##----------------
  }
}

###################################################################
#####plot
# par(mfrow=c(1,2),mai=c(0.8,0.8,1,0.1))
plot(apply(resultsS,1,mean), ylim=c(0,1),xaxt="n",main = bquote("Unpaired test, n="~.(na)),ylab = expression("Power"),xlab = expression(delta),type="b",pch=19)
axis(1,at=seq(1,ldeltas,by=1),labels = deltas)
lines(apply(resultsrankS,1,mean), pch=3,col="red")
lines(apply(resultsF,1,mean), pch=5,col="blue")
lines(apply(resultsRoger,1,mean), pch=9,col="purple")
abline(c(alpha, 0), lty=2)
legend("topleft",legend = c(expression(paste(infinity,"-S/median test")),expression(paste(infinity,"-rankS test")), "F-test",expression(paste(chi^2,"-test"))),col=c("black","red","blue","purple"),lty=c(1,1,1,1),pch=c(19,NA,NA,NA))#pch = c(1,3,5,9)

# save(precalculated,resultsS,resultsrankS,resultsF,resultsRoger,Sout_MCH0,rankSout_MCH0,file=paste("n",na,"unpaired.RData",sep = ""))

