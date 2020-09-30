#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//' @description
//' description
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
long long int chooseRec(int n, int k) {
  if (k == 0 || k == n)
    return 1;
  return chooseRec(n - 1, k - 1) + chooseRec(n - 1, k);
}

// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//' @description
//' description
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
long long int choose(int n, int k)
{
  long long int res = 1;

  // Since C(n, k) = C(n, n-k)
  if(k > n - k){
    k = n - k;
  }

  // Calculate value of
  // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (int i = 0; i < k; ++i){
    res *= (n - i);
    res /= (i + 1);
  }

  return res;
}
