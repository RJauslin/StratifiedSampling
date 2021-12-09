


rm(list = ls())


library(sampling)
data("swissmunicipalities")
eps <- 1e-7 # epsilon tolerance
n <- 200 # sample size
pik <- inclusionprobabilities(swissmunicipalities$POPTOT,n)
mask <- (pik < (1 - eps)) & (pik > eps)
pik <-  pik[mask]
# pik <- pik[sample(1:length(pik))]
# UPmaxentropy(pik)


plot(density(pik))

var(rgamma(1000,1,10))


plot(sort(rgamma(100,1,5)))





rm(list = ls())
eps <- 1e-7 # epsilon tolerance
set.seed(1)
pik <- inclusionprobabilities(rgamma(1000,0.2,10),100)
mask <- (pik < (1 - eps)) & (pik > eps)
pik <-  pik[mask]



psipik2 <- function(n,pik){
  w <- pik/(1-pik)
  
  out = rep(0,length(pik))
  for(i in 1:n){
    tmp = w*(1.0 - out)
    if(any( (1.0 - out) > 1 - 1e-16 ) && i > 1 ){ # not 1-out that goes out of the scope
      break;
    }
    if(any( (i * tmp/sum(tmp)) > 1 - 1e-16 ) && i > 1 ){ # not 1-out that goes out of the scope
      print("FALSE")
      break;
    }
    out = i * tmp/sum(tmp)
    
  }
  return(out)
}




test <- psipik2(30,pik)
max(test)


n <- sum(pik)
for(nn in 1:n){
  cat("Step: ",nn," min :",min(psipik2(nn,pik)) , "max :",max(psipik2(nn,pik)),"\n\n")
  if(min(psipik2(nn,pik)) < 1e-16 || max(psipik2(nn,pik)) > (1-1e-16)){
    break;
  }
}






# psirec <- function(n,pik){
#   if(n == 1){
#     return(rep(0,length(pik)))
#   }else{
#     w <- pik/(1-pik)
#     tmp <- w*(1.0 - psirec(n-1,pik))
#     return( n* tmp/sum(tmp))
#   }
# }


# psirec(45,pik)

# logit(pik)
# 
# min(pik/(1-pik))
# max(pik/(1-pik))
# 
# psi(n,pik/(1-pik))


pik <- pik[-(pik < 1e-2)]
psipik2 <- function(n,pik){
  w <- pik/(1-pik)
  which(w < 1e-3)
  
  out = rep(0,length(pik))
  for(i in 1:n){

    # print(min(1.0-out))
    # print(max(1.0-out))
    tmp = w*(1.0 - out)
    # Sys.sleep(1)
    # plot(tmp)
    out = i * tmp/sum(tmp)
    if(any(out > (1 - 1e-8))){
      
    }
  }
  return(out)
}

min(psipik2(1,pik))


test <- psipik2(61,piktmp)




sum(test) 
min(test)
max(test)
w <- pik/(1-pik)
tmp <- w*(1.0 - test)

max(62*tmp/sum(tmp))
min(62*tmp/sum(tmp))




tmp < (1/61)*sum(tmp)

pik <- piktmp
for(nn in 1:n){
  cat("Step: ",nn," min :",min(psipik2(nn,pik)) , "max :",max(psipik2(nn,pik)),"\n\n")
  if(min(psipik2(nn,pik)) < 1e-8 || max(psipik2(nn,pik)) > (1-1e-8)){
    break;
  }
}

sum(UPpoisson(pik))
