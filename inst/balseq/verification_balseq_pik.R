

N=100
n <- 10
p=10
# pik=runif(N)
pik=inclusionprobabilities(runif(N),n)
Xaux=array(rnorm(N*p,3,1),c(N,p))

Xspread <- cbind(runif(N),runif(N))
Xaux=cbind(pik,Xaux)
# Xaux <- cbind(pik)

s <- balseq(pik,Xspread,Xaux)


SIM=1000
pikh=rep(0,N)
for(i in 1:SIM){
  print(i)
  pikh=pikh+balseq(pik,Xspread,Xaux)
} 
pikh=pikh/SIM

pik
pikh

t=(pikh-pik)/sqrt(pik*(1-pik)/SIM)
sum(abs(t)<1.96)/N


