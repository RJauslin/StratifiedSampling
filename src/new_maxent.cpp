#include <Rcpp.h>
#include <random>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector psi(int n,NumericVector& w) {
  NumericVector out = w/sum(w);
  for(int i = 2;i <= n; i++){
    out = i* w*( 1.0 - out )/sum( w* (1.0 - out) );
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector psipik(int n,NumericVector& pik) {
  NumericVector w = pik/(1-pik);
  NumericVector out(pik.size());
  
  for(int i = 1;i <= n; i++){
    // Rcout << out << std::endl;
    out = i* w*( 1.0 - out )/sum( w* (1.0 - out) );
  }
  return(out);
}


// [[Rcpp::export]]
NumericVector logit(NumericVector x){
  return(log(x/(1-x)));
}



// [[Rcpp::export]]
NumericVector w(NumericVector& pik,double tol = 1e-6,int max_iter = 500){
  int n = round(sum(pik));
  NumericVector w_tmp = pik/(1.0-pik);
  NumericVector pikstar = psi(n,w_tmp);
  double err = sum(abs(pik - pikstar));
  int counter = 1;
  while(err > tol && counter < max_iter){
    Rcout << min(w_tmp) << std::endl;
    w_tmp = w_tmp* (psi(n,w_tmp)/(1-psi(n,w_tmp))) * (1-pik)/pik;
    counter++;
  }
  return(w_tmp);
}



// [[Rcpp::export]]
NumericVector piktilde(NumericVector& pik, double tol = 1e-6,int max_iter = 500){
  int n = round(sum(pik));
  NumericVector pikstar = clone(pik);
  NumericVector tmp = pikstar/(1.0-pikstar);
  double err = sum(abs(pikstar - psi(n,tmp)));// [[Rcpp::export]]
  
  int  counter = 1;
  while(err > tol && counter < max_iter ){
    Rcout << min(pikstar) << std::endl;
    pikstar = pikstar + (pik - psi(n,tmp));
    tmp = pikstar/(1.0-pikstar);
    err = sum(abs(pik - psi(n,tmp)));
    counter++;
 } 
  return(pikstar);
  
}


// [[Rcpp::export]]
NumericVector lambda(NumericVector piktilde){
  NumericVector out = log(piktilde/(1-piktilde));
  return(out);
}




// [[Rcpp::export]]
Rcpp::IntegerVector sample_int(int n,int N) {
  Rcpp::IntegerVector pool = Rcpp::seq(1, N);
  std::random_shuffle(pool.begin(), pool.end());
  return pool[Rcpp::Range(0, n-1)];
} 



/*** R



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

logit(pik)

min(pik/(1-pik))
max(pik/(1-pik))

psi(n,pik/(1-pik))

for(nn in 1:n){
  cat("Step: ",nn," min :",min(psipik(nn,pik)) , "max :",max(psipik(nn,pik)),"\n\n")
  if(min(psipik(nn,pik)) < 1e-8 || max(psipik(nn,pik)) > (1-1e-8)){
    break;
  }
}

w <- w_(pik,sum(pik),length(pik))



# x <- seq(0,1,0.001)
# plot(x,log(x/(1-x)))


plot(pik,log(tmp))

min(psi(sum(pik),pik/(1-pik)))
max(psi(sum(pik),pik/(1-pik)))
pikt <- piktilde(pik,max_iter = 100)
w <- pikt/(1.0 - pikt)
qfromw(w,sum(pik))



library(sampling)
# set.seed(1)
N <- 100
n <- 10
pik <- inclusionprobabilities(runif(N),n)
sum(pik)



psi(n,pik/(1-pik))
psipik(n,pik)


min(psi(sum(pik),pik/(1-pik)))
max(psi(sum(pik),pik/(1-pik)))

tmp <- pik/(1-pik)

plot(pik,tmp)

pikt <- piktilde(pik)
w <- pikt/(1.0 - pikt)
qfromw(w,sum(pik))


# phi(n,pik/(1-pik))
# phi2(n,pik/(1-pik))


system.time(piktilde(pik))
system.time(UPMEpiktildefrompik(pik))

test <- piktilde(pik)
l <- lambda(test)



sum(UPMEpiktildefrompik(pik))


*/
