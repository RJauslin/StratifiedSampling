rm(list = ls())
N <- 1000

x1 <- rgamma(N,4,25)
x2 <- rgamma(N,4,25)
x3 <- rgamma(N,4,25)
x4 <- rgamma(N,4,25)

nstrata <- 25

strata <- as.matrix(rep((1:nstrata),each = (N/nstrata))) 
Xcat <- disjMatrix(strata)

#-------- CASE 1 pik equal and integer in each strata

pik <- inclusionprobastrata(strata,rep(2,nstrata)) #↓ one by strata
sum(pik)
X <- as.matrix(matrix(c(x1,x2),ncol = 2))

# 
# s <- fbs(X,strata,pik)

 y1 <- 20*strata + rnorm(1000,120)
 y2 <- 500 + 5*x1 + 5*x2 + rnorm(1000,0,270)
 y3 <- 500 + 100*x1 + 100*x2 + 100*x3 + 100*x4 + rnorm(1000,0,1000)
 y4 <- 500 + 200*x1 + 100*x2 + 100*x3 + 50*x4 + rnorm(1000,0,1000)
 
 
 SIM = 1000
 t1 <- matrix(rep(0,2*SIM),ncol = 2)
 varApp(X,strata,pik,y1)
 
 
 
 for(i in 1:SIM){
   print(i)
   
    
   # s <- stratifiedcube(cbind(pik,X),strata,pik)
   # s <- fbs(X,strata,pik)
   s <- balstrat(X,strata,pik)
   y <- y1
   y_ht <- sum(y[which(s==1)]/pik[which(s == 1)])
   t1[i,1] <- (sum(y_ht) - sum(y))^2
   t1[i,2] <- varEst(X,strata,pik,s,y)
 }
 
 
colSums(t1)/SIM
varApp(X,strata,pik,y1)

