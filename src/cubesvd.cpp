#include <RcppArmadillo.h>
#include "projOp.h"

// arma::vec onestep_cpp(arma::vec u,arma::vec pik,double EPS=0.0000001){
//   
//   arma::uword N = pik.size();
//   double l1 = 1e+200;
//   double l2 = 1e+200;
//   double l = 1e-9;
//   
//   for(arma::uword k = 0; k < N; k++){
//     if(u[k]> 0){
//       l1 = std::min(l1,(1.0 - pik[k])/u[k]);
//       l2 = std::min(l2,pik[k]/u[k]);
//     }
//     if(u[k]< 0){
//       l1 = std::min(l1,-pik[k]/u[k]);
//       l2 = std::min(l2,(pik[k]-1.0)/u[k]);
//     }
//   }
//   if(Rcpp::runif(1)[0]<l2/(l1+l2)){
//     l = l1;
//   }else{
//     l = -l2;
//   }
//   for(arma::uword k = 0; k < N; k++){
//     pik[k] = pik[k] + l*u[k];
//     if(pik[k] < EPS){
//       pik[k] = 0;
//     }
//     if(pik[k] > (1-EPS)){
//       pik[k] = 1;
//     }
//   }
//   return(pik);
// }


void onestep_cpp(arma::vec& u,arma::vec& pik,double EPS=0.0000001){
  
  arma::uword N = pik.size();
  double l1 = 1e+200;
  double l2 = 1e+200;
  double l = 1e-9;
  
  for(arma::uword k = 0; k < N; k++){
    if(u[k]> 0){
      l1 = std::min(l1,(1.0 - pik[k])/u[k]);
      l2 = std::min(l2,pik[k]/u[k]);
    }
    if(u[k]< 0){
      l1 = std::min(l1,-pik[k]/u[k]);
      l2 = std::min(l2,(pik[k]-1.0)/u[k]);
    }
  }
  if(Rcpp::runif(1)[0]<l2/(l1+l2)){
    l = l1;
  }else{
    l = -l2;
  }
  for(arma::uword k = 0; k < N; k++){
    pik[k] = pik[k] + l*u[k];
    if(pik[k] < EPS){
      pik[k] = 0;
    }
    if(pik[k] > (1-EPS)){
      pik[k] = 1;
    }
  }
  // return(pik);
}



// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//'
//'
//' @param X matrix of auxiliary variables.
//' @param pik vector of inclusion probabilities
//' @param EPS tolerance
//'
//' @return a sample
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::vec ffphase_cpp(arma::mat& X,arma::vec pik,double EPS=0.0000001){

  arma::mat D = arma::diagmat(1/pik);
  arma::mat A = D*X;
  unsigned int J = X.n_cols; // initial size of B
  
  arma::uvec i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first"); // find first index of B
  arma::mat B = (A.rows(i)).t(); // extract B of A
  
  while(i.size() > 0){
    // std::cout << i.size() << std::endl;
    
    arma::mat kern = arma::null(B);
    arma::uword N = pik.size();
    arma::vec u(N);
    u = kern.col(0);
    
    arma::vec pik_tmp = pik.elem(i);
    // pik.elem(i) =  onestep_cpp(u,pik.elem(i));
    onestep_cpp(u,pik_tmp);
    pik.elem(i) = pik_tmp;
    i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first");
    B = (A.rows(i)).t();
    // std::cout << B.n_cols << std::endl;
    // std::cout << B.n_rows << std::endl << std::endl;
    if(i.size() < (J+1)){
      arma::mat kern = arma::null(B);
      if(kern.empty()){
        break;
      }
    }
  }
  
  return(pik);
}


// 
// 
// // [[Rcpp::depends(RcppArmadillo)]]
// //' @title title
// //'
// //'
// //'
// //' @param X matrix of auxiliary variables.
// //' @param pik vector of inclusion probabilities
// //' @param EPS tolerance
// //'
// //' @return a sample
// //'
// //' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
// //'
// //' @export
// // [[Rcpp::export]]
// arma::vec landing_svd(arma::mat X,arma::vec pikstar,double EPS=0.0000001){
//   
//   arma::mat D = arma::diagmat(1/pikstar);
//   arma::mat A = D*X;
//   unsigned int J = X.n_cols; // initial size of B
//   
//   arma::uvec i = arma::find(pikstar > EPS && pikstar < (1-EPS), J+1, "first"); // find first index of B
//   arma::mat B = (A.rows(i)).t(); // extract B of A
// 
//   arma::mat U;
//   arma::vec s;
//   arma::mat V;
//   
//   while(i.size() > 1){
//     arma::mat one = arma::ones<arma::mat>(i.size(),1);
//     
//     // SVD 
//     arma::svd_econ(U,s,V,B,"right","dc");
//     arma::vec u(i.size());
//     u = V.col(V.n_cols - 1);
//     u = u - projOp(u,one);
//     pikstar.elem(i) =  onestep_cpp(u,pikstar.elem(i));
//     i = arma::find(pikstar > EPS && pikstar < (1-EPS));
//     B = (A.rows(i)).t();
//   }
//   
//   // LAST UNIT CHOOSE BY USING BERNOULLI VARIABLE
//   if(i.size() == 1){
//     pikstar[i[0]] = R::rbinom(1,pikstar[i[0]]);
//   }
//   
//   return(pikstar);
// }
// 
// 
// 
// // [[Rcpp::depends(RcppArmadillo)]]
// //' @title title
// //'
// //'
// //'
// //' @param X matrix of auxiliary variables.
// //' @param pik vector of inclusion probabilities
// //' @param EPS tolerance
// //'
// //' @return a sample
// //'
// //' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
// //'
// //' @export
// // [[Rcpp::export]]
// arma::vec cubesvd(arma::mat X,arma::vec pik,double EPS=0.0000001){
//   
//   arma::mat D = arma::diagmat(1/pik);
//   arma::mat A = D*X;
//   unsigned int J = X.n_cols; // initial size of B
//   
//   arma::uvec i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first"); // find first index of B
//   arma::mat B = (A.rows(i)).t(); // extract B of A
//   
//   while(i.size() > 0){
//     // std::cout << i.size() << std::endl;
//     
//     arma::mat kern = arma::null(B);
//     arma::uword N = pik.size();
//     arma::vec u(N);
//     u = kern.col(0);
//     
//     pik.elem(i) =  onestep_cpp(u,pik.elem(i));
//     i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first");
//     B = (A.rows(i)).t();
//     // std::cout << B.n_cols << std::endl;
//     // std::cout << B.n_rows << std::endl << std::endl;
//     if(i.size() < (J+1)){
//       arma::mat kern = arma::null(B);
//       if(kern.empty()){
//         break;
//       }
//     }
//   }
//   
//   // LANDING PHASE BY USING SMALLEST VECTOR
//   arma::mat U;
//   arma::vec s;
//   arma::mat V;
//   
//   while(i.size() > 1){
//     arma::mat one = arma::ones<arma::mat>(i.size(),1);
// 
//     // SVD 
//     arma::svd_econ(U,s,V,B,"right","dc");
//     arma::vec u(i.size());
//     u = V.col(V.n_cols - 1);
//     u = u - projOp(u,one);
//     pik.elem(i) =  onestep_cpp(u,pik.elem(i));
//     i = arma::find(pik > EPS && pik < (1-EPS));
//     B = (A.rows(i)).t();
//   }
// 
//   // LAST UNIT CHOOSE BY USING BERNOULLI VARIABLE
//   if(i.size() == 1){
//     pik[i[0]] = R::rbinom(1,pik[i[0]]);
//   }
// 
//   return(pik);
// }
// 


/*** R
library(sampling)
rm(list = ls())
N = 3000
n = 20
p = 20
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))



library("microbenchmark")

micro_timings = microbenchmark(ffphase(X,pik),
                               ffphase_cpp(X,pik),
                               BalancedSampling::flightphase(pik,X),
                               times = 30)
plot(micro_timings)
micro_timings






rm(list = ls())
set.seed(1)
eps <- 1e-13
library(Matrix)
N <- 250
Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),70),
                sampling::inclusionprobabilities(runif(N),50),
                sampling::inclusionprobabilities(runif(N),30)),ncol = 3)
X <- PM(Pik)$PM
pik <- PM(Pik)$P
dim(X)
order = 2
method = 1
comment = TRUE
EPS = 1e-11
system.time(test <- ffphase_cpp(X,pik))
# t(X)%*%test
# t(X)%*%pik
system.time(pikstar <- fastflightcubeSPOT(X, pik))
# t(X)%*%pikstar
# t(X)%*%pik
system.time(pikstar2 <- sampling::fastflightcube(X, pik))
# t(X)%*%pikstar2
# t(X)%*%pik
system.time(test <- BalancedSampling::flightphase(pik,X))
system.time(test <- fast.flight.cube(X,pik))
t(X)%*%test
t(X)%*%pik

*/
