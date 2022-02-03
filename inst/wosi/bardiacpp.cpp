#include <RcppArmadillo.h>
#include "inclprob.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


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
  arma::vec pika(pikr.begin(),N,false);
  
  // deep copy
  arma::vec pik(pika);
  
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
  
  arma::uvec s_tmp = osod(pik_tmp);
  m = pikb;
  
    
  // // MAIN LOOP
  // int c = 0;
  // for(int i = k;i < N; i++){
  //   c = c+1;
  //   m = m + pik[i];
  //   
  //   // no need to do a step if equal to 0 or 1
  //   if(m >= H){
  //     
  //     double pikai = H - (m - pik[i]);
  //     arma::vec pikstar(c);
  //     for(int l = 0;l < c-1; l++){
  //       pik_tmp[l] = pik[l];
  //     }
  //     pik_tmp[k] = pikai;
  //     arma::vec inclprob(pikstar,H);
  //     
  // 
  //     
  //     // index = arma::regspace<arma::uvec>(i+1,i + bound);
  // 
  //     
  //    c = 0; 
  //   }
  //   
  //   c++;
  // }
  
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