#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Binomial Coefficient
//'
//' @description
//' It calculates the number of sets with size \eqn{k} that can be chosen from a set of size \eqn{n}.
//'
//' @param n An integer. (greater than k)
//' @param k An integer.
//'
//' @return An integer "n choose k".
//'
//' @author Raphael Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @seealso
//' \code{\link[base:choose]{base::choose}}
//'
//' @examples
//' choose(10,5)
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
