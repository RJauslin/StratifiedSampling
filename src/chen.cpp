#include <Rcpp.h>
#include <random>
using namespace Rcpp;



// [[Rcpp::export]]
double T(NumericVector& w,int i,int j){
  NumericVector tmp = pow(w,i);
  double res = pow(w[j-1],i);
  return(sum(tmp) - res);
}

// [[Rcpp::export]]
double RR(NumericVector& w,int k,int j){
  NumericVector T_tmp(k);
  for(int i = 1;i <= k;i++){ 
    T_tmp[i-1] = T(w,i,j);
  }
  double beta = 0.0;
  NumericVector R_tmp(k);
  for(int i = 1; i <= k;i++){
    for(int l = 1;l <= i;l++){
      // beta = (1.0/i)*pow(-1.0,l+1)*T(w,l,j);
      beta = (1.0/i)*pow(-1.0,l+1)*T_tmp[l-1];
      if((i-l)==0){
        R_tmp[i-1] = R_tmp[i-1] + beta;
      }else{
        R_tmp[i-1] = R_tmp[i-1] + beta * R_tmp[i-l-1];  
      }
    }
    // Rcout << R_tmp << std::endl;
  }
  return(R_tmp[k-1]);
}


// [[Rcpp::export]]
NumericVector w_(NumericVector& pik,int n, int J,int max_iter = 50, double eps = 1e-6){
  int N = pik.size();
  NumericVector w = clone(pik);
  // NumericVector w(N);
  // std::copy( w.begin(), w.end(), pik.begin() );
  
  double num(0.0);
  double den(0.0);
  double cond(0.0);
  NumericVector w_tmp(N);
  int t = 1;
  do{
    
    for(int kk = 0; kk < N; kk++){
      w_tmp[kk] = w[kk];
    }
    for(int i = 1;i <= J; i++){
      num = RR(w,n-1,J);
      den = RR(w,n-1,i);
      w[i-1] = pik[i-1]*(num/den);
    }
    cond =  max(abs(w/w_tmp - 1.0));
    Rcout << cond << std::endl;
  } while (cond > eps && t < max_iter);
  
  return(w);
}



/*** R
library(sampling)
rm(list = ls())
pik <- inclusionprobabilities(seq(1,100,by = 1),10)
N <- length(pik)
n <- sum(pik)
w <- pik
 # w <- w_(pik,n,N)

# verification
RR(pik,sum(pik),length(pik))
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
# pik <- inclusionprobabilities(runif(100),10)
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

s/SIM
pik


####################################################33

rm(list = ls())
set.seed(1)
library(sampling)
N <- 2000
n <- 30
pik <- inclusionprobabilities(runif(N),n)
any(pik > (1-1e-7))
system.time(w <- w_(pik,n,N))
q <- qfromw(w,n)
s <- sfromq(q)



rm(list = ls())
library(sampling)
data("swissmunicipalities")
eps <- 1e-7 # epsilon tolerance
n <- 200 # sample size
pik <- inclusionprobabilities(swissmunicipalities$POPTOT,n)
mask <- (pik < (1 - eps)) & (pik > eps)
pik <-  pik[mask]

n <- sum(pik)
N <- length(pik)
system.time(w <- w_(pik,n,N))






*/
