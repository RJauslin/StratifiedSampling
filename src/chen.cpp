#include <Rcpp.h>
#include <random>
using namespace Rcpp;



// [[Rcpp::export]]
double T(NumericVector w,int i,int j){
  NumericVector tmp = pow(w,i);
  double res = pow(w[j-1],i);
  return(sum(tmp) - res);
}

// [[Rcpp::export]]
double RR(NumericVector w,int k,int j){
  NumericVector T_tmp(k);
  for(int i = 1;i <= k;i++){
    T_tmp[i-1] = T(w,i,j);
  }
  double beta = 0.0;
  NumericVector R_tmp(k);
  for(int i = 1; i <= k;i++){
    for(int l = 1;l <= i;l++){
      beta = (1.0/i)*pow(-1.0,l+1)*T(w,l,j);
      if((i-l)==0){
        R_tmp[i-1] = R_tmp[i-1] + beta;
      }else{
        R_tmp[i-1] = R_tmp[i-1] + beta * R_tmp[i-l-1];  
      }
    }
  }
  return(R_tmp[k-1]);
}

// [[Rcpp::export]]
NumericVector w_(NumericVector pik,int n, int N){
  NumericVector w = clone(pik);
  double num(0.0);
  double den(0.0);
  for(int t = 1; t <= 10; t++){
    for(int i = 1;i <= N; i++){
      num = RR(w,n-1,N);
      den = RR(w,n-1,i);
      w[i-1] = pik[i-1]*(num/den);
    }
  }
  return(w);
}



/*** R
library(sampling)

pik <- inclusionprobabilities(seq(1,100,by = 1),10)
N <- length(pik)
n <- sum(pik)

 w <- w_(pik,n,N)

# verification
RR(w,3,N-1)
T(w,3,N-1)

R(w,3,N-1)
fT(w,3,N-1)
 


0.5*fT(w,1,N-1)^2 - 0.5*fT(w,2,N-1) # second entry
(1/6)*fT(w,1,N-1)^3 - (1/6)*fT(w,1,N-1)*fT(w,2,N-1) - (1/3)*fT(w,2,N-1)*fT(w,1,N-1) + (1/3)*fT(w,3,N-1) # thrid entry

pik <- c(0.1,0.4,0.7,0.8)
w_(pik,2,4)


updatew(pik,n,N)


pik <- inclusionprobabilities(seq(1,100,by = 1),10)




############ TEST


rm(list = ls())
pik <- inclusionprobabilities(seq(1,100,by = 1),10)
pik <- inclusionprobabilities(runif(100),10)
N <- length(pik)
n <- sum(pik)
#w <- updatew(pik,n,N)
w <- w_(pik,n,N)
# w <- w - mean(w)
eps <- 1e-8
q <- qfromw(w,n)
s <- rep(0,length(pik))
SIM <- 200000
for(i in 1 :SIM){
  tmp <- sfromq(q)  
  if(any(abs(sum(tmp) - sum(pik)) > eps)){
    cat("error")
    break;
  }
  s <- s + tmp
}

p <- s/SIM
pik 

s[which((p-pik)/sqrt(pik*(1-pik)/SIM) > 2)]/SIM
pik[which((p-pik)/sqrt(pik*(1-pik)/SIM) > 2)]
sum(s)







*/
