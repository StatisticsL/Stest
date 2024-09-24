# install_github("StatisticsL/Stest/Stest")
library(Stest)
library(quantmod)
AMZN_data <- getSymbols(Symbols = "AMZN", src = "yahoo", from = "2022-08-25", to = "2024-08-27", auto.assign = FALSE)

timeseries=diff(log(AMZN_data$AMZN.Close))
y=as.vector(timeseries[-1])

n=length(y)
p=n
A="diff"
b=rep(0,n-1)
X=diag(n)
alpha=0.05
tau=0.1
M=10000

#####------------
precalculated=precalculatedmatrices(A=A,n=n)
out_rescaled=rescaleAb(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,precalculated=precalculated)
homopowerrescaling=out_rescaled$homopowerrescaling
Sout_MCH0=MCH0(X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,rescale=T,precalculated = precalculated,homopowerrescaling = homopowerrescaling,q=1)

obj=Stest(y=y,X=X,A=A,b=b,tau=tau,alpha=alpha,M=M,Rank=F,rescale = T,outMCH0=Sout_MCH0, precalculated = precalculated,homopowerrescaling = homopowerrescaling,q=1,nlambda = 10)
plot(obj,main =expression("Last 2 years AMZN's daily log-returns"),xlab = expression("Days"),ylab = expression("Log-returns"))

save(precalculated,out_rescaled,homopowerrescaling,Sout_MCH0,obj,file="financialtimeseries.RData")


