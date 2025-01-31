#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Disjunctive
//' @name disj
//' @description
//' This function transforms a categorical vector into a matrix of indicators.
//'
//' @param strata_input A vector of integers that represents the categories.
//'
//' @return A matrix of indicators.
//'
//'
//' @examples
//' strata <- rep(c(1,2,3),each = 4)
//' disj(strata)
//'
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix disj(Rcpp::IntegerVector strata_input) {
  
  arma::uvec strata = Rcpp::as<arma::uvec>(strata_input);
  int N = strata.size();
  arma::uvec cat = arma::unique(strata);
  int ncat = cat.size();

  arma::uvec val = arma::regspace<arma::uvec>(0,1,ncat-1); // 0 to unique(strata)
  arma::uvec strata_tmp = strata;
  for(int i = 0;i<ncat;i++){
    strata_tmp.replace(cat[i],val[i]);
  }

  
  arma::umat m(N,ncat,arma::fill::zeros);
  for(int i = 0;i < N;i++){
    m(i,strata_tmp(i)) = 1;
  }
  // return(m);
  
  Rcpp::IntegerMatrix out = Rcpp::wrap(m);
  return(out);
}


/***R

strata <- rep(c(1,2,3),each = 4)
disj(strata)


strata=c(-2,3,-2,3,4,4,4,-2,-2,3,4,0,0,0)

strata <- rep(c(1,2,3),each = 4)
disj(strata)



set.seed(1)
strata <- sample(x = 1:6, size = 50, replace = TRUE)
system.time(M <- disj(strata))
system.time(M <- model.matrix(~as.factor(strata)-1))


N <- 100000
strata <- sample(x = 1:400, size = N, replace = TRUE)
system.time(M <- disj(strata))
system.time(M <- model.matrix(~as.factor(strata)-1))
system.time(M <- sampling::disjunctive(strata))


*/


// [[Rcpp::depends(RcppArmadillo)]]
//' @title Number of categories
//' @name ncat
//' @description
//' This function returns the number of factor in each column of a categorical matrix.
//'
//' @param Xcat_input A matrix of integers that contains categorical vector in each column.
//'
//' @return A row vector that contains the number of categories in each column.
//'
//'
//' @export
//' 
//' @examples
//' Xcat <-  matrix(c(sample(x = 1:6, size = 100, replace = TRUE),
//'             sample(x = 1:6, size = 100, replace = TRUE),
//'             sample(x = 1:6, size = 100, replace = TRUE)),ncol = 3)
//' ncat(Xcat)
// [[Rcpp::export]]
Rcpp::NumericVector ncat(Rcpp::IntegerMatrix Xcat_input){
  arma::umat Xcat = Rcpp::as<arma::umat>(Xcat_input);
  int p = Xcat.n_cols;
  arma::rowvec m(p,arma::fill::zeros);
  for(int i = 0; i < p; i++){
    arma::uvec cat = unique(Xcat.col(i));
    int tmp = cat.size();
    m(i) = tmp;
  }
  
  Rcpp::NumericVector out = Rcpp::wrap(m);
  // out.attr("dim") = R_NilValue;
  return(out);
  
}

/***R

strata <-  sample(x = 1:6, size = N, replace = TRUE)
ncat(strata)


Xcat <-  matrix(c(sample(x = 1:6, size = 100, replace = TRUE),
                sample(x = 1:6, size = 100, replace = TRUE),
                sample(x = 1:6, size = 100, replace = TRUE)),ncol = 3)
ncat(Xcat)

for(i in 1:100){
  Xcat <- cbind(Xcat, sample(x = 1:6, size = 100, replace = TRUE))
}
system.time(n <- apply(Xcat,MARGIN = 2,FUN <- function(x){nlevels(as.factor(x))}))
system.time(test <- ncat(Xcat))

*/



arma::vec ncat_internal(arma::umat Xcat){
  int p = Xcat.n_cols;
  arma::vec out(p,arma::fill::zeros);
  for(int i = 0; i < p; i++){
    arma::uvec cat = unique(Xcat.col(i));
    int tmp = cat.size();
    out(i) = tmp;
  }
  return(out);
}

arma::umat disj_internal(arma::uvec strata) {
  
  int N = strata.size();
  arma::uvec cat = arma::unique(strata);
  int ncat = cat.size();
  
  arma::uvec val = arma::regspace<arma::uvec>(0,1,ncat-1); // 0 to unique(strata)
  arma::uvec strata_tmp = strata;
  for(int i = 0;i<ncat;i++){
    strata_tmp.replace(cat[i],val[i]);
  }
  
  arma::umat m(N,ncat,arma::fill::zeros);
  for(int i = 0;i < N;i++){
    m(i,strata_tmp(i)) = 1;
  }
  return(m);
}


// [[Rcpp::depends(RcppArmadillo)]]
//' @title Disjunctive for matrix  
//' @name disjMatrix
//' @description
//' This function transforms a categorical matrix into a matrix of indicators variables.
//'
//' @param strata_input A matrix of integers that contains categorical vector in each column.
//'
//' @return A matrix of indicators.
//'
//' 
//' @examples
//' Xcat <-  matrix(c(sample(x = 1:6, size = 100, replace = TRUE),
//'             sample(x = 1:6, size = 100, replace = TRUE),
//'             sample(x = 1:6, size = 100, replace = TRUE)),ncol = 3)
//' disjMatrix(Xcat)
//'
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix disjMatrix(Rcpp::IntegerMatrix strata_input) {

  
  arma::umat strata = Rcpp::as<arma::umat>(strata_input);
  
  
  arma::uword N = strata.n_rows;
  arma::uword p = strata.n_cols;
  arma::vec all_cat = ncat_internal(strata);
  
  int n_all_cat = 0;
  for(arma::uword i = 0;i < all_cat.size(); i++){
    n_all_cat += all_cat[i];
  }
  
  arma::umat m(N,n_all_cat,arma::fill::zeros);
  arma::vec subind = arma::round(cumsum(all_cat)-1);
  arma::uvec ind = arma::conv_to<arma::uvec>::from(subind);

  for(arma::uword i = 0;i < p;i++){
    arma::uvec tmp = strata.col(i);
    arma::umat tmp_mat(N,all_cat(i),arma::fill::zeros);
    tmp_mat = disj_internal(tmp);
    
    if(i == 0){
      m.cols(0, ind(i)) = tmp_mat;
    }else{
      m.cols(ind(i-1)+1, ind(i)) = tmp_mat;
    }

  }

  // return(m);
  
  Rcpp::IntegerMatrix out = Rcpp::wrap(m);
  // out.attr("dim") = R_NilValue;
  return(out);
  
}


/***R

Xcat <-  matrix(c(sample(x = 1:6, size = 100, replace = TRUE),
             sample(x = 1:6, size = 100, replace = TRUE),
             sample(x = 1:6, size = 100, replace = TRUE)),ncol = 3)
disjMatrix(Xcat)


N <- 1000
Xcat <- sample(x = 1:40, size = N, replace = TRUE)
for(i in 1:1000){
  Xcat <- cbind(Xcat,sample(x = 1:40, size = N, replace = TRUE))
}
# sum(ncat(Xcat))
# ncat(Xcat)


system.time(test <- apply(as.matrix(Xcat),MARGIN = 2,FUN <- function(x){as.matrix(model.matrix(~as.factor(x)-1))}))
system.time(test <-  disjMatrix(Xcat))


disjMatrix(as.matrix(Xcat[,1]))

*/


