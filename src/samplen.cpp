#include <RcppArmadillo.h>
#include <algorithm>
#include <string>
#include <iostream>
#include <vector>
#include <functional>
#include "choose.h"


// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//' @description
//' description Rosetta code combn
//'
//'
//' @param x x
//'
//' @details
//'
//' details
//'
//' @return a vector
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @seealso
//' func
//'
//' @examples
//'
//' @export
// [[Rcpp::export]]
arma::mat samplen(int N, int n){

  arma::uword count = choose(N,n);
  arma::mat out(N,count,arma::fill::zeros);
  arma::uvec tmp(n,arma::fill::zeros);
  arma::vec tmp2(N,arma::fill::zeros);

  std::string bitmask(n, 1); // K leading 1's
  bitmask.resize(N, 0); // N-K trailing 0's

  // print integers and permute bitmask
  int c = 0;
  do {
    int c1 = 0;
    for(int i = 0; i < N; ++i){
      if(bitmask[i]){
        tmp(c1) = i;
        c1++;
      }
    }
    tmp2(tmp) += 1;
    out.col(c) = tmp2;
    tmp2(tmp) -= 1;
    c++;
  } while(std::prev_permutation(bitmask.begin(), bitmask.end()));
  return(out);
}

/*** R

system.time(test <- IneqCube::choose(10,3))
system.time(test <- base::choose(10,3))


system.time(test <- samplen(26,10))
system.time(test2 <- combn(26,10))
system.time(test <- sampling::writesample(10,26))
test
test2

*/
