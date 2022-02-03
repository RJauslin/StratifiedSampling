


N <- 1000
n <- 100
pik <- inclusionprobabilities(runif(N),n)


t <- c()
for(i in 1:10000){
  s <- UPpoisson(pik)
  t[i] <- sum(s)
}

sum(t)/10000


rm(list = ls())
N <- 1000
n <- 100
pik <- inclusionprobabilities(runif(N),n)


SIM <- 1000

t <- c()
tmp <- rep(0,N)
for(i in 1:SIM){
  k <- 0
  s <- c()
  print(i)
  while(k <= N){
    # print(k)
    if(k+1 == N){
      print(pik[N])
      s[N] <- round(pik[N])
      break;
    }
    b <- c_bound(pik[(k+1):N])
    
    pik_tmp <- pik[(k+1):(k+b)]
    # print(sum(pik_tmp))
    s_tmp <- osod(pik_tmp)
    # print(s_tmp)
    s <- c(s,s_tmp)
    k <- k + b
  }
  print(length(s))
  tmp <- tmp + s
  t <- c(t,sum(s))
}
tmp <- tmp/SIM
sum(t)/SIM
test <- (tmp-pik)/sqrt(pik*(1-pik)/SIM)
sum(abs(test)<1.96)/N
