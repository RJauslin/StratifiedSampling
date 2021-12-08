#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector phi(int n,NumericVector& w) {
  if(n == 1){
    return w/sum(w);
  }else{
    return n * w * (1.0 - phi(n - 1,w))/(sum(w*(1.0 - phi(n - 1,w))));
  }
}

// [[Rcpp::export]]
NumericVector piktilde(NumericVector& pik, double tol = 1e-6,int max_iter = 500){
  int n = round(sum(pik));
  Rcout << n << std::endl;
  NumericVector pikstar = clone(pik);
  NumericVector tmp = pikstar/(1.0-pikstar);
  double err = sum(abs(pikstar - phi(n,tmp)));// [[Rcpp::export]]
  
  int  counter = 1;
  while(err > tol && counter < max_iter ){
    Rcout << err << std::endl;
    pikstar = pikstar + (pik - phi(n,tmp));
    tmp = pikstar/(1.0-pikstar);
    err = sum(abs(pik - phi(n,tmp)));
    counter++;
 } 
  return(pikstar);
  
}







/*** R

library(sampling)
set.seed(1)
N <- 500
n <- 16
pik <- inclusionprobabilities(runif(N),n)
sum(pik)
# system.time(phi(n,pik/(1-pik)))

system.time(piktilde(pik))
system.time(UPMEpiktildefrompik(pik))

sum(UPMEpiktildefrompik(pik))


*/
