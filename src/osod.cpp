#include <RcppArmadillo.h>
#include "inclprob.h"
#include "c_bound.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


//' @title One-step One Decision sampling method
//'
//' @description This function implements the One-step One Decision method. It can be used using equal or unequal inclusion probabilities. The method is particularly useful for selecting a sample from a stream. 
//' 
//' @param pikr A vector of inclusion probabilities.
//' 
//' @details
//' 
//' The method sequentially transforms the vector of inclusion probabilities into a sample whose values are equal to 0 or 1. The method respects the inclusion probabilities and can
//' handle equal or unequal inclusion probabilities.
//' 
//' The method does not take into account the whole vector of inclusion probabilities by having a sequential implementation. This means that the method is fast and can be implemented in a flow.
//' 
//' @return A vector with elements equal to 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for rejected units.
//'
//' @author Raphael Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @seealso \code{\link{c_bound}}
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
IntegerVector osod(NumericVector pikr){ 
  
  
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
      
      bound = c_bound(pik.subvec(i,N-1));
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

N <- 1000
n <- 100
pik <- inclprob(runif(N),n)
s <- osod(pik)

*/
