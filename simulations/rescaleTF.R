# install_github("StatisticsL/Stest/Stest")
library(Stest)
###----------------------------------------------
n=100
p=20
m=2
# X=matrix(rt(n*p, df=2),n,p)/sqrt(n)
X=matrix(rnorm(n*p),n,p)
beta=rep(0,p)
A=diag(m)
A=cbind(A,matrix(0,m,p-m))
A[1,]=A[1,]*3
b=rep(0,m)

alpha=0.05
tau=0.8
M=10000
#######################################################
precalculated=precalculatedmatrices(A=A,n=n)
out_rescaled=rescaleAb(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,precalculated=precalculated)
homopowerrescaling=out_rescaled$homopowerrescaling
Sout_MCH0_T=MCH0(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,rescale=T,precalculated = precalculated,homopowerrescaling = homopowerrescaling,q=1)
Sout_MCH0_F=MCH0(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,rescale=F,precalculated = precalculated,q=1)



numiter=50
deltas=seq(0,1.6,length=8); ldeltas=length(deltas)
resultsS_T1=matrix(NA,ldeltas,numiter)
resultsS_F1=matrix(NA,ldeltas,numiter)


for(j in 1:ldeltas){
  print(j)
  for (iter in 1:numiter) {
    err=rnorm(n)
    delta=c(deltas[j],rep(0,p-1))
    # delta=c(0,deltas[j],rep(0,p-2))
    y=X%*%(beta+delta)+err

    resultsS_T1[j,iter]=Stest(y=y,X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,rescale=T,outMCH0=Sout_MCH0_T,homopowerrescaling = homopowerrescaling,precalculated = precalculated,q=1,nlambda = 0)$test_value
    resultsS_F1[j,iter]=Stest(y=y,X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,rescale=F,outMCH0=Sout_MCH0_F, precalculated = precalculated,q=1,nlambda = 0)$test_value

  }
}


resultsS_T=matrix(NA,ldeltas,numiter)
resultsS_F=matrix(NA,ldeltas,numiter)

for(j in 1:ldeltas){
  print(j)
  for (iter in 1:numiter) {
    err=rnorm(n)
    # delta=c(deltas[j],rep(0,p-1))
    delta=c(0,deltas[j],rep(0,p-2))
    y=X%*%(beta+delta)+err

    resultsS_T[j,iter]=Stest(y=y,X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,rescale=T,outMCH0=Sout_MCH0_T,homopowerrescaling = homopowerrescaling, precalculated = precalculated,q=1,nlambda = 0)$test_value
    resultsS_F[j,iter]=Stest(y=y,X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,rescale=F,outMCH0=Sout_MCH0_F, precalculated = precalculated,q=1,nlambda = 0)$test_value

  }
}


Wout_rescaled=rescaleW(X=X,A=A,b=b,homopowerrescaling = homopowerrescaling,tau=tau,alpha=alpha,M=M,Rank=F,precalculated=precalculated)

####plot----------------------------------------------------
par(mfrow=c(2,2),mai=c(0.7,0.8,0.4,0.8))
boxplot(Wout_rescaled$allwvec,main=expression("Not rescaled"),xaxt="n")
axis(1,at=c(1,2),labels = c(expression(W[1]),expression(W[2])))
boxplot(Wout_rescaled$allwvecrescaled,main=expression("Rescaled"),xaxt="n")
axis(1,at=c(1,2),labels = c(expression(W[1]),expression(W[2])))

plot(apply(resultsS_T1,1,mean),ylab="Power", xlab=expression(delta),ylim=c(0,1),xaxt="n",main = expression(H[paste(1,",",1,sep="")]),type="b",pch=19)
axis(1,at=seq(1,ldeltas,by=1),labels = round(deltas,2))
lines(apply(resultsS_F1,1,mean),pch=3,col="red")
abline(c(alpha, 0), lty=2)


plot(apply(resultsS_T,1,mean),ylab="Power", xlab=expression(delta),ylim=c(0,1),xaxt="n",main =  expression(H[paste(1,",",2,sep="")]),type="b",pch=19)
axis(1,at=seq(1,ldeltas,by=1),labels = round(deltas,2))
lines(apply(resultsS_F,1,mean),pch=3,col="red")
abline(c(alpha, 0), lty=2)
legend(x=4,y=0.4,legend = c(expression("Rescaled"), expression("Not rescaled")),col=c("black","red"),lty=c(1,1),pch = c(19,NA))

# save(precalculated,homopowerrescaling,Sout_MCH0_T,Sout_MCH0_F,resultsS_T1,resultsS_F1,resultsS_T,resultsS_F,
#      Wout_rescaled,file = "rescaleTF.Rdata")
