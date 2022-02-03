#include <RcppArmadillo.h>
#include "c_bound.h"
#include "inclprob.h"
#include "osod.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


arma::uvec osod_internal(arma::vec pik,
                   bool full = false){ 
  
  
  // transfer memory for arma vector
  double eps = 1e-8;
  int N = pik.size();
  
  // sum of inclusion probabilities (double and int to check if it is equal to integer or not)
  double n = arma::sum(pik);
  double n_int = std::round(n);
  
  
  // ghost unit, depending if the sum of the inclusion probabilities are equal to integer a ghost unit is added.
  bool ghost = false;
  if(abs(n-n_int) > 1e-6){
    pik.resize(N+1);
    N = pik.size();
    pik[N-1] = ceil(n)-n;
    n = arma::sum(pik);
    ghost = true;
  }
  
  
  // initialization of temporary variables
  int bound = 0;
  double n_tmp = 0.0;
  arma::uvec index;
  arma::vec pik_s;
  double w = 0.0;
  
  // MAIN LOOP
  for(int i = 0;i < N-1; i++){
    
    w = pik[i]; // weight
    
    // no need to do a step if equal to 0 or 1
    if(w > eps && w < (1.0-eps)){
      
      if(full == true){
        bound = N-1-i;
      }else{
        bound = c_bound(pik.subvec(i,N-1));  
      }
      
      index = arma::regspace<arma::uvec>(i+1,i + bound);
      
      pik_s = pik.elem(index);
      n_tmp = arma::sum(pik_s) + w;
      
      
      pik[i] = 0.0;// put unit i equal to 0
      pik_s = inclprob(pik_s,n_tmp); // update inclusion prob
      
      
      if(runif(1)[0] < w){
        
        pik.elem(index) = (pik.elem(index) - pik_s*(1.0-w))/w;
        pik[i] = 1.0;
      }else{
        pik.elem(index) = pik_s;
      }
      
    }
  }
  
  // rounding and return as IntegerVector
  if(ghost == true){
    arma::uvec s(N-1);
    for(int i = 0;i< N-1;i++){
      s[i] = pik[i];
    }
    return(s);
  }else{
    arma::uvec s(N);  
    for(int i = 0;i< N;i++){
      if(pik[i] > (1-eps)){
        s[i] = 1;
      }
    }
    return(s);
  }
  
}











//' @title wosi
//'
//' @description sakdhfa
//' 
//' @param pikr A vector of inclusion probabilities.
//' 
//' @details
//' 
//' 
//' @return A vector with elements equal to 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for rejected units.
//'
//' @author Raphael Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' 
//' @examples
//' 
//' N <- 1000
//' n <- 100
//' pik <- inclprob(runif(N),n)
//' s <- osod(pik)
//' 
//' @export
// [[Rcpp::export]]
IntegerVector wosicpp(NumericVector pikr,
                   int H){ 
  
  
  // transfer memory for arma vector
  double eps = 1e-8;
  int N = pikr.size();
  arma::vec pikarma(pikr.begin(),N,false);
  
  // deep copy
  arma::vec pik(pikarma);
  
  // sum of inclusion probabilities (double and int to check if it is equal to integer or not)
  double n = arma::sum(pik);
  double n_int = std::round(n);
  
  
  // ghost unit, depending if the sum of the inclusion probabilities are equal to integer a ghost unit is added.
  bool ghost = false;
  if(abs(n-n_int) > 1e-6){
    pik.resize(N+1);
    N = pik.size();
    pik[N-1] = ceil(n)-n;
    n = arma::sum(pik);
    ghost = true;
  }
  
  
  arma::uvec out(N);
  

  // FIRST STEP
  int k = 0;
  double m = pik[k];
  do{
    k++;
    m += pik[k];
  } while (m < H);
  
  
  double pika = H - (m  - pik[k]);
  double pikb = pik[k] - pika;

  arma::vec pik_tmp(k);
  for(int l = 0;l < (k-1);l++){
    pik_tmp[l] = pik[l];
  }
  pik_tmp[k] = pika;

  arma::uvec s_tmp = osod_internal(pik_tmp);
  m = pikb;
  Rcout << s_tmp << std::endl;
    
  // MAIN LOOP
  int c = 0;
  for(int i = k;i < N; i++){

    m = m + pik[i];
    c++;
    // no need to do a step if equal to 0 or 1
    if(m >= H){

      double pikai = H - (m - pik[i]);
      arma::vec pikstar(c);
      for(int l = 0;l < c-1; l++){
        pikstar[l] = pik[l];
      }
      pikstar[c] = pikai;
      arma::vec pikstar2 = inclprob(pikstar,H);
      
      Rcout << *(s_tmp.end()-1) << std::endl;


      // index = arma::regspace<arma::uvec>(i+1,i + bound);


     c = 0;
    }
  }

  // rounding and return as IntegerVector
  if(ghost == true){
    IntegerVector s(N-1);
    for(int i = 0;i< N-1;i++){
      s[i] = pik[i];
    }
    return(s);
  }else{
    IntegerVector s(N);  
    for(int i = 0;i< N;i++){
      if(pik[i] > (1-eps)){
        s[i] = 1;
      }
    }
    return(s);
  }
  
}




/*** R


rm(list = ls())
library(sampling)
pik<-c(.4,.7,.1,.2,.4,.2,.2,.2,.1,.6,.2,.1,.4,.5,.5,.6,.1,.2)
pik <- inclusionprobabilities(runif(100),10)
H <- 2
s <- wosicpp(pik,H)
s

*/


