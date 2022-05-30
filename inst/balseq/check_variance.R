rm(list = ls())
N <- 1000
n <- 30

Xaux <- matrix(rgamma(N*5,4,25),ncol = 5)
Xspread <- as.matrix(cbind(runif(N),runif(N)))
pik <- rep(n/N,N)
y <- Xaux%*%c(1,1,1,1,1) + rnorm(N,0,1)



SIM <- 10000

v_estimated <- v2 <- v3 <- v4 <- rep(0,N)
for(i in 1:SIM){
  print(i)
  s <- BalancedSampling::lcube(pik,Xspread,Xaux)
  s_srs <- which(sampling::srswor(n,N) == 1)
  s_01 <- rep(0,N)
  s_01[s] <- 1
  
  v_estimated[i] <- vEst(Xaux,pik,y,s_01)
  v2[i] <- varB(Xaux[s,],Xspread[s,],pik[s],y[s])
  v3[i] <- varDBS(Xaux[s,],Xspread[s,],pik[s],y[s])
  v4[i] <- varDBS(cbind(pik[s]),Xspread[s,],pik[s],y[s])
}
mean(v_estimated)
mean(v2)
mean(v3)
mean(v4)
