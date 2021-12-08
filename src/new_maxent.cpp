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
NumericVector piktilde(NumericVector& pik, double tol = 1e-6,int max_iter = 500){
  int n = round(sum(pik));
  NumericVector pikstar = clone(pik);
  NumericVector tmp = pikstar/(1.0-pikstar);
  double err = sum(abs(pikstar - psi(n,tmp)));// [[Rcpp::export]]
  
  int  counter = 1;
  while(err > tol && counter < max_iter ){
    Rcout << pikstar << std::endl;
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
pik <- pik[sample(1:length(pik))]
# UPmaxentropy(pik)

pikt <- piktilde(pik)
w <- pikt/(1.0 - pikt)
qfromw(w,sum(pik))



library(sampling)
set.seed(1)
N <- 100
n <- 10
pik <- inclusionprobabilities(runif(N),n)
sum(pik)

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
