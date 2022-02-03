

rm(list = ls())

library(sampling)
 pik<-c(.3,.4,.5,.2,.1,.3,.3,.4,.3,.6,.2,.1,.2,.5,.4,.6,.5,.1)
# pik<-c(.4,.7,.1,.2,.4,.2,.2,.2,.1,.6,.2,.1,.4,.5,.5,.6,.1,.2)


# pik<-inclusionprobabilities(runif(100),10)
pik
cumsum(pik)
H <- 2



SIM=10000
N<-length(pik)
pikh=rep(0,N)
PI=array(0,c(N,N))
for(i in 1:SIM)
{
  a1=wosi2(pik,H)
  a1
  print(sum(a1))
  pikh=pikh+a1
  PI=PI+a1%*%t(a1)
}
pikh=pikh/SIM
PI=PI/SIM
t=(pikh-pik)/sqrt(pik*(1-pik)/SIM)

print(sum(abs(t)<1.96)/N)
