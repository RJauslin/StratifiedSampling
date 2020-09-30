#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
//' @title disjunctive vector
//'
//' @description
//' This function transform a categorical variable into a matrix of indicators.
//'
//' @param strata A vector of integer that represents the category.
//'
//' @return A matrix of indicators.
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @examples
//'   strata <- rep(c(1,2,3),each = 4)
//'   disj(strata)
//'
//' @export
// [[Rcpp::export]]
arma::umat disj(arma::uvec strata) {
  int N = strata.size();
  arma::uvec cat = arma::unique(strata);
  int ncat = cat.size();

  arma::uvec val = arma::regspace<arma::uvec>(0,1,ncat-1); // 0 to unique(strata)
  arma::uvec strata_tmp = strata;
  for(arma::uword i = 0;i<ncat;i++){
    strata_tmp.replace(cat[i],val[i]);
  }

  arma::umat m(N,ncat,arma::fill::zeros);
  for(arma::uword i = 0;i < N;i++){
    m(i,strata_tmp(i)) = 1;
  }
  return(m);
}


/***R
strata=c(-2,3,-2,3,4,4,4,-2,-2,3,4,0,0,0)

strata <- rep(c(1,2,3),each = 4)
disj(strata)


rm(list = ls())
set.seed(1)
strata <- sample(x = 1:6, size = 50, replace = TRUE)
system.time(M <- disj(strata))
system.time(M <- model.matrix(~as.factor(strata)-1))


rm(list = ls())
N <- 100000
strata <- sample(x = 1:400, size = N, replace = TRUE)
system.time(M <- disj(strata))
system.time(M <- model.matrix(~as.factor(strata)-1))
system.time(M <- sampling::disjunctive(strata))


*/



// [[Rcpp::depends(RcppArmadillo)]]
//' @title number of category
//'
//' @description
//' This function returns the number of factor in each column of a categorical matrix.
//'
//' @param Xcat A matrix of integer that contains categorical variable in column.
//'
//' @return A row vector that contains the number of category in each column.
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::rowvec ncat(arma::umat Xcat){
  int p = Xcat.n_cols;
  arma::rowvec out(p,arma::fill::zeros);
  for(arma::uword i = 0; i < p; i++){
    arma::uvec cat = unique(Xcat.col(i));
    int tmp = cat.size();
    out(i) = tmp;
  }
  return(out);
}

/***R
rm(list = ls())
Xcat <-  sample(x = 1:6, size = 100, replace = TRUE)
for(i in 1:10000){
  Xcat <- cbind(Xcat, sample(x = 1:6, size = 100, replace = TRUE))
}
system.time(n <- apply(Xcat,MARGIN = 2,FUN <- function(x){nlevels(as.factor(x))}))
system.time(test <- ncat(Xcat))

*/



// [[Rcpp::depends(RcppArmadillo)]]
//' @title disjunctive Matrix
//'
//' @description
//' This function transform a categorical matrix into a matrix of indicators variables.
//'
//' @param A matrix of integer that contains categorical variable in column.
//'
//' @return A matrix of indicators.
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::umat disjMatrix(arma::umat strata) {

  int N = strata.n_rows;
  int p = strata.n_cols;
  arma::rowvec all_cat = ncat(strata);
  int n_all_cat = sum(all_cat);

  arma::umat m(N,n_all_cat,arma::fill::zeros);
  arma::rowvec subind = arma::round(cumsum(all_cat)-1);
  arma::uvec ind = arma::conv_to<arma::uvec>::from(subind);


  for(arma::uword i = 0;i < p;i++){
    arma::uvec tmp = strata.col(i);
    arma::umat tmp_mat(N,all_cat(i),arma::fill::zeros);
    tmp_mat = disj(tmp);

    if(i == 0){
      m.cols(0, ind(i)) = tmp_mat;
    }else{
      m.cols(ind(i-1)+1, ind(i)) = tmp_mat;
    }

  }

  return(m);
}


/***R
rm(list = ls())

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


// [[Rcpp::depends(RcppArmadillo)]]
//' @title findBarma
//'
//' @description
//' findB
//'
//' @param X a Matrix
//' @param Xcat a Matrix
//'
//' @return a matrix
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::mat findBarma(arma::mat X,
                    arma::umat Xcat){
  // double eps = 1e-9;
  int pInit = X.n_cols;
  int N = X.n_rows;

  arma::umat Xcat_tmp = Xcat.rows(0,pInit-1);
  int n_all_cat = sum(ncat(Xcat_tmp));
  int n_all_cat_tmp = 0;
  int p = 0;


  while(n_all_cat != n_all_cat_tmp){
    n_all_cat_tmp = n_all_cat;
    p = pInit + n_all_cat;
    if(p > N){
      p = N;
    }
    // std::cout << p << std::endl;
    Xcat_tmp = Xcat.rows(0,p-1);

    n_all_cat = sum(ncat(Xcat_tmp));
  }

  arma::umat Xdev = disjMatrix(Xcat_tmp);
  arma::mat Xdev_tmp = arma::conv_to<arma::mat>::from(Xdev);
  arma::mat final = arma::join_rows(X.rows(0,p-1),Xdev_tmp);
  return(final);

}



/***R
rm(list = ls())

library(sampling)
N <- 5000
Xcat <-as.matrix(data.frame(cat1 = rep(1:40,each = N/40),
                  cat2 = rep(1:50,each = N/50),
                  cat3 = rep(1:500,times = 2)))
                  # cat3 = rep(1:100,each = N/100)))

pik <- inclusionprobastrata(Xcat[,1],rep(1,40))
p <- 30
X <- as.matrix(matrix(rnorm(N*p),ncol = 30))
dim(findB(X,Xcat))

system.time(B1 <- findBarma(X,Xcat))
system.time(B2 <- findB(X,Xcat))


rm(list = ls())

A <- matrix(c(-0.7105,   4.1029,  -1.8246,
              -2.0019 , -0.5236 ,  4.1916,
              4.3996,  -2.8083  ,-0.5213,
              0.6630 , -2.8091  , 7.5276,
              -2.9156 , -0.7766 , -1.1729,
              -0.3270  , 4.2548 ,  5.3200,
              0.1400  ,-4.0253  , 2.8248,
              -4.2845  , 2.7870 ,  0.7671,
              -3.0131 , -0.2912 ,  7.2370,
              0.6044  , 4.3184  , 3.1805,
              -1.9504 ,  3.2579 ,  0.2756,
              -1.6467 ,  5.7019 , -6.5195,
              1.4168  ,-1.0080  ,-4.7708,
              1.0052  ,-5.1903  , 2.9107),ncol = 3,byrow = T)

Xcat <- matrix(c(4,         5  ,       9,
                 4,         5  ,       9,
                 4,         4  ,       8,
                 3 ,        4  ,       8,
                 4 ,        5  ,      10,
                 2 ,        2  ,       3,
                 4 ,        5  ,      10,
                 1 ,        2  ,       3,
                 4 ,        5  ,      10,
                 4 ,        5  ,      10,
                 4 ,        5  ,      10,
                 4 ,        4  ,       8,
                 4 ,        5 ,       10,
                 4 ,        5,        10),ncol = 3,byrow = T)


findBarma(A,Xcat)

A <- matrix(c(-0.0804,  -6.1586,   2.8135,
1.0696,  -0.7477,   3.6091,
-1.9891,   3.7850,  -4.4040,
0.7696,   0.9360,   3.0155,
0.1584,  -2.2411,   2.7247,
1.6843,   2.3359,  -4.9651,
5.1468,  -2.0724,  -2.7368,
0.9597,  -0.0283,   2.4393,
1.2555,  -7.7195,  -5.2076,
6.4425,   0.5949,  -6.3745,
4.0282,   3.0615,  -7.6527,
-3.6085,  -1.3441,  -0.0273,
-4.1969,  -4.0039,   3.8278),ncol = 3, byrow = T)


Xcat <- matrix(c(1,         2,         3,
3,         4 ,        8,
2 ,        2  ,       3,
4 ,        5  ,      10,
3 ,        4  ,       7,
4 ,        5,        10,
4 ,        5  ,      10,
4 ,        5  ,      10,
4 ,        5  ,      10,
4 ,        5  ,      10,
3 ,        4  ,       7,
4 ,        5  ,      10,
4 ,        5  ,      10),ncol = 3,byrow = T)


findBarma(A,Xcat)

*/
