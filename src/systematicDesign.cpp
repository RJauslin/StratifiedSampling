#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
//' @encoding UTF-8
//' @title systematic desgin
//'
//' @description
//'
//' Find all possible systematic sample from a vector of inclusion probabilities \code{pik}.
//'
//' @param pik vector of inclusion probabilities..
//'
//' @return A matrix of size at most N x N with elements equal to 0 or 1. The value 1 indicates that the uni is selected while the value 0 is for non-chosen unit.
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' @export
// [[Rcpp::export]]
Rcpp::List systematicDesign(arma::vec pik){

  //Value to round value like -1e-17...
  double eps = 1e-13;
  double R = 1e+9;




  //index of pik equal 1 or 0
  arma::uvec ones = arma::find(pik > (1-eps));
  arma::uvec zeros = arma::find(pik < eps);

  // add difference between ceiling(pik) and sum(pik) (possibly 0)
  arma::vec pik1(1);
  pik1.fill(ceil(arma::sum(pik)) - sum(pik));
  arma::vec piks = join_cols(pik, pik1);

  // Number of element
  int N = piks.n_elem;

  // vk and r
  arma::vec vk = arma::cumsum(piks);


  //round value otherwise -1e-17 can appears
  vk = round(R*vk)/R;
  arma::vec vk1 = vk - arma::floor(vk);
  arma::vec r = arma::sort(vk1);

  // add 1 to r
  arma::vec one(1);
  one.fill(1);
  r = join_cols(r, one);

  // centered value and p (rounded p to be sure that is correctly equal to 0)
  arma::vec cent = (r(arma::span(0, N-1)) + r(arma::span(1, N)))/2;
  arma::vec p = r(arma::span(1, N)) - r(arma::span(0, N-1));
  p = round(R*p)/R;

  // add 0 to vk
  arma::vec zero(1);
  zero.fill(0);
  vk = join_cols(zero,vk);

  // loop that select sample
  arma::vec A(N+1);
  arma::umat final(N,N);
  arma::vec tmp3(N+1);
  for(int i = 0;i < N;i++){
    tmp3.fill(cent(i));
    tmp3 = round(R*(vk-tmp3))/R;
    A = tmp3 - arma::floor(tmp3);
    final.col(i) = A(arma::span(0, N-1)) >  A(arma::span(1, N));
  }

  // remove empty last line
  final.shed_row(N-1);

  //remove duplicated value
  arma::uvec dupl = find(p < eps);
  final.shed_cols(dupl);
  final = final.t();
  arma::uvec dupl_inv = find(p >= eps);
  p = p.elem(dupl_inv);

  // be sure that pik with 1 or 0 equal 1 or 0 on each column
  arma::uvec tmp4(final.n_rows,arma::fill::ones);
  for(unsigned int k = 0;k< ones.size();k++){
    final.col(ones(k)) = tmp4;
  }
  for(unsigned int j = 0;j < zeros.size();j++){
    tmp4.fill(0.0);
    final.col(zeros(j)) = tmp4;
  }


  // return(final);
  return Rcpp::List::create(Rcpp::Named("probas")=p,
                            Rcpp::Named("samples") = final);
}



/*** R

# pik <- runif(10000)
# system.time(test1 <- systematicDesign(pik))
# system.time(test2 <- systematic.find.design(pik))


rm(list = ls())
library(sampling)
# pik <- inclusionprobabilities(runif(10),8)

# pik with 1 and 0
pik <- c(0.3, 1.0, 0.0, 0.7, 0.4, 0.15, 0.25, 0.2)
sum(pik)
# same without 0 and 1
# pik <- c(0.3, 0.7, 0.4, 0.15, 0.25, 0.2)
# sum(pik)

pik <- rep(0.3,12)
sum(pik)

system.time(test1 <- systematicDesign(pik))
test1
rowSums(test1$samples)
system.time(test2 <- SystematicDesign(pik))
test2
rowSums(test2)


round(test1,9) == round(test2,9)



pik <- rep(1,10)
systematicDesign(pik)
SystematicDesign(pik)


*/
