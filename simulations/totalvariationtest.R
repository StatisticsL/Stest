# install_github("StatisticsL/Stest/Stest")
library(Stest)
##############################################################################################################################
n=10
A="diff"
b=rep(0,n-1)
X=diag(n)
alpha=0.05
tau=0.8
M=10
#######################################################
precalculated=precalculatedmatrices(A=A,n=n)
Sout_MCH0=MCH0(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,precalculated = precalculated,q=1)
rankSout_MCH0=MCH0(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=T,precalculated = precalculated,q=1)

###############################################################################################################################
numiter=10
deltas=seq(0,2,by=0.2)
ldeltas=length(deltas)
###########################################################
resultsS=matrix(NA,ldeltas,numiter)
resultsrankS=matrix(NA,ldeltas,numiter)
df=10
for (j in 1:ldeltas) {
  print(j)
  for (iter in 1:numiter) {
    y1=rt(n/2,df=df)*sqrt((df-2)/df)
    y2=rt(n/2,df=df)*sqrt((df-2)/df)*(1+deltas[j])
    y=c(y1,y2)

    resultsS[j,iter]=Stest(y=y,X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,outMCH0=Sout_MCH0, precalculated = precalculated,q=1,nlambda = 0)$test_value
    resultsrankS[j,iter]=Stest(y=y,X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=T,outMCH0=rankSout_MCH0,precalculated = precalculated,q=1,nlambda = 0)$test_value
  }
}

###############################################################################################################################
###############################################################################################################################

###left plot----------------------------------------------
par(mfrow=c(1,2),mai=c(0.8,0.8,1,0.1))
plot(apply(resultsS,1,mean),ylab="Power", xlab=expression(delta),ylim=c(0,1),xaxt="n",main = expression("Student errors, df=10"),type="b",pch=19)
axis(1,at=seq(1,ldeltas,by=1),labels = deltas)
lines(apply(resultsrankS,1,mean),pch=3,col="red")
abline(c(alpha, 0), lty=2)
legend("topleft",legend = c(expression(paste(infinity,"-S test")),  expression(paste(infinity,"-rankS test"))),col=c("black","red"),lty=c(1,1),pch = c(19,NA))
####################################################################
# save(precalculated,resultsS,resultsrankS,Sout_MCH0,rankSout_MCH0,file="totalvariationtest.RData")

####################################################################
###right plot----------------------------------------------
y1=rt(n/2,df=df)*sqrt((df-2)/df)
y2=rt(n/2,df=df)*sqrt((df-2)/df)*(1+2)
y=c(y1,y2)

aa=Stest(y=y,X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,outMCH0=Sout_MCH0, precalculated = precalculated,q=1)$test_value
if(aa==1){
  lambda0data=ztfaql(y=y,X=X,A=A,b=b,tau=tau,Rank=F,outtype="2",precalculated=precalculated,q=1)
  plot(y,ylab=expression(y[i]),main=expression(paste("Example for rejected ",H[0])),xlab=expression(i))
  lines(rep(lambda0data$indexmax+0.5,2),range(y),lty=2)
}
###############################################################################################################################

