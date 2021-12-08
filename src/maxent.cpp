#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
NumericMatrix qfromw(NumericVector& wr,const int& n){

  arma::vec w(wr.begin(),wr.size(),false);
  
  int N = w.size();
  arma::mat expa(N,n,arma::fill::zeros);
  for(int i = 0;i < N;i++){
    expa(i,0) = arma::sum(w.subvec(i,N-1));
  }
  for(int i = 1;i < n; i++){
    expa(N-i-1,i) = exp(arma::sum(log(w.subvec(N-i-1,N-1))));
  }
  for(int i = N-3; i >= 0; i--){
    for(int z = 1; z < std::min(N-i,n); z++){
      expa(i,z) = w[i]*expa(i+1,z-1) + expa(i+1,z);
    }
  }
  
  NumericMatrix q(N,n);
  for(int i = 0;i < N;i++){
    q(i,0) = w[i]/expa(i,0);
  }
  for(int i = 1;i < n; i++){
    q(N-i-1,i) = 1;
  }
  for(int i = N-3; i >= 0; i--){
    for(int z = 1; z < std::min(N-i,n); z++){
      q(i,z) = w[i]*expa(i+1,z-1)/expa(i,z);
    }
  }
  return(q);
}

//' @export
// [[Rcpp::export]]
IntegerVector sfromq(const NumericMatrix& q){
  int N = q.nrow();
  int n = q.ncol();
  IntegerVector s(N);
 
  for(int k = 0;k < N;k++){
    if(n != 0){
      if(runif(1)[0] < q(k,n-1)){
        s[k] = 1;
        n = n-1;
      }
    }
  }
  return(s);
}

//' @export
// [[Rcpp::export]]
NumericVector pikfromq(NumericMatrix& qr){
  int N = qr.nrow();
  int n = qr.ncol();
  
  arma::mat q(qr.begin(),N,n,false);
  arma::mat pro(N,n,arma::fill::zeros);
  
  pro(0,n-1) = 1.0;
  for(int i = 1;i< N;i++){
    for(int j = 1;j< n;j++){
      pro(i,j) += pro(i-1,j)*(1-q(i-1,j));
      pro(i,j-1) +=  pro(i-1,j)*(q(i-1,j));
    }
  }
  for(int i = 1;i < N; i++){
    pro(i,0) += pro(i-1,0)*(1-q(i-1,0));
  }
  arma::vec out = arma::sum(pro%q,1);
  NumericVector out2(N);
  for(int i = 0; i < N;i++){
    out2[i] = out[i];
  }
  return(out2);
}

//' @export
// [[Rcpp::export]]
NumericVector piktfrompik(NumericVector& pik){

  
 int N = pik.size();
 int n = round(sum(pik));
 NumericVector pikt(Rcpp::clone(pik));
 double arr = 1.0;
 double eps = 1e-8;

 NumericVector w(N);
 NumericMatrix q(N,n);
 NumericVector pikt1(N);

 while(arr > eps){
  w = pikt/(1.0 - pikt);
  q = qfromw(w,n);
  pikt1 = pikt + pik - pikfromq(q);
  arr = sum(abs(pikt - pikt1));
  pikt = pikt1;
 }
 return(pikt);

}


//' @title Maximum entropy sampling
//'
//' @description Maximum entropy sampling with fixed sample size. It can handle unequal inclusion probabilities.
//' 
//' @param pikr A vector of inclusion probabilities.
//' 
//' @details
//' The sampling design maximizes the entropy design:
//' \deqn{I(p) = - \sum s p(s) log[p(s)].}
//' 
//' This function is a C++ implementation of \code{\link[sampling:UPmaxentropy]{UPmaxentropy}}.
//' More details could be find in Tille (2006).
//' @return A vector with elements equal to 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for rejected units.
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' 
//' @references 
//' Tille, Y. (2006), Sampling Algorithms, springer
//' 
//' @export
// [[Rcpp::export]]
IntegerVector maxent(NumericVector& pikr){
  
  double eps = 1e-6;
  int N = pikr.size();
  arma::vec pik_tmp(pikr.begin(),pikr.size(),false);
  
  // find not equl to 0 value
  arma::uvec i = arma::find(pik_tmp < 1-eps);
  
  // arma to Numeric vector 
  arma::vec pik_tmp2 = pik_tmp.elem(i);
  int N_tmp = pik_tmp2.size();
  NumericVector pik(N_tmp);
  for(int j = 0;j< N_tmp;j++){
    pik[j] = pik_tmp2[j];
  }

  // all step computation
  int n = round(sum(pik));
  NumericVector pikt = piktfrompik(pik);
  NumericVector w = pikt/(1-pikt);
  NumericMatrix q = qfromw(w,n);
  IntegerVector s2 = sfromq(q);
  
  
  //put right 0-1 value at the right place
  IntegerVector s(N);
  s.fill(1);
  for(int j = 0; j < N_tmp;j++){
    s[i[j]] = s2[j];
  }
  
  return(s);
}

//' @export
// [[Rcpp::export]]
NumericMatrix pik2frompik(NumericVector pikr, NumericVector wr){
  
  
  double eps = 1e-6;
  int N = pikr.size();
  arma::vec pik(pikr.begin(),pikr.size(),false);
  arma::vec w(wr.begin(),wr.size(),false);
  int n = arma::sum(pik);
  
  NumericMatrix M(N,N);
  
  
  for(int k = 0; k < N; k++){
    for(int l = 0; l < N; l++){
      if( (abs(pik[k] - pik[l]) > eps) & (k != l)){
      // if(pik[k] != pik[l] & k != l){
        M(k,l) = (pik[k]* w[l] - pik[l]*w[k])/(w[l] - w[k]);
      }else{
        M(k,l) = -1.0;
      }
    }
  }
  
  for(int i = 0; i < N; i++){
    M(i,i) = pik[i];
  }
  
  
  double tt = 0.0;
  int comp = 0;
  double cc = 0.0;
  for(int k = 0; k < N; k++){
    tt = 0.0;
    comp = 0;
    for(int l = 0; l < N; l++){
      // if(M(k,l) != -1.0){
      if(abs(M(k,l) + 1.0) > eps){  
        tt = tt + M(k,l);  
      }else{
        comp = comp + 1.0;
      }
    }
    cc = (n * pik[k] - tt)/comp;
    for(int l = 0; l < N; l++){
      if(abs(M(k,l) + 1) < eps){
      // if(M(k,l) == -1.0){
        M(k,l) = cc;
      }
    }  
  }
  
  return(M);
  
  
}

/*** R
rm(list = ls())
library(sampling)
data("swissmunicipalities")
eps <- 1e-7 # epsilon tolerance
n <- 200 # sample size
pik <- inclusionprobabilities(swissmunicipalities$POPTOT,n)
mask <- (pik < (1 - eps)) & (pik > eps)
pik <-  pik[mask]
pik <- pik[sample(1:length(pik))]
# UPmaxentropy(pik)

v <- function(M){
  require(Matrix)
  image(as(M,"sparseMatrix"))
}


w <- pik/(1-pik)


# pik <- c(0.07,0.17,0.41,0.61,0.83,0.91)
# UPMEpiktildefrompik(pik)
# pikt <- c(0.1021,0.2238,0.4417,0.5796,0.7794,0.8734)
# pik/(1-pik)
# pikt/(1-pikt)
# w <- c(0.116,0.295,0.810,1.411,3.614,7.059)
# n <
# Uqf <- function (w, n) 
# {
  N = length(w)
  expa = array(0, c(N, n))
  # fill first column
  for (i in 1:N){
    expa[i, 1] = sum(w[i:N])
  } 
  min(expa[,1]) # all filled
  
  # fill diagonal bottom left to right up
  for (i in (N - n + 1):N){
    cat(i,N-i + 1,"\n\n")
    expa[i, N - i + 1] = exp(sum(log(w[i:N])))
    print(expa[i, N - i + 1])
  }
  v(expa)
  
  for (i in (N - 2):1){
    for (z in 2:min(N - i, n)) {
      Sys.sleep(1)
      # cat(i,z," use ",i+1,z-1,"and ",i+1,z,"\n\n")
      expa[i, z] = (w[i]*expa[i + 1, z - 1]) + expa[i + 1, z]
      print(expa[i,z])
      # if(expa[i,z] < 1e-20){
        # print(expa[i,z])  
      # }
      
    }
  } 
  v(expa)
  
  q = array(0, c(N, n))
  for (i in N:1) q[i, 1] = w[i]/expa[i, 1]
  for (i in N:(N - n + 1)) q[i, N - i + 1] = 1
  v(q)
  
  for (i in (N - 2):1){
    for (z in 2:min(N - i, n)){
      q[i, z] = w[i] * expa[i + 1, z - 1]/expa[i, z]
    } 
  } 
  
  print(any(is.na(q)))
  q
}


test2 <- Uqf(w,sum(pik))
any(which(is.na(test2)))


test <- qfromw(pik/(1-pik),sum(pik))


piktfrompik(pik)[length(pik)]
maxent(pik)






























pik=c(0.07,0.17,0.41,0.61,0.83,0.91)

n=sum(pik)
pikt=UPMEpiktildefrompik(pik)
pikt=piktfrompik(pik)
w=pikt/(1-pikt)
q=UPMEqfromw(w,n)
UPMEsfromq(q)

pik2frompik(pik,w)
UPMEpik2frompikw(pik,w)

*/




//' @title Joint inclusion probabilities of maximum entropy.
//'
//' @description This function computes the matrix of the joint inclusion of the maximum entropy sampling with fixed sample size. It can handle unequal inclusion probabilities.
//' 
//' @param pikr A vector of inclusion probabilities.
//' 
//' @details
//' The sampling design maximizes the entropy design:
//' \deqn{I(p) = - \sum s p(s) log[p(s)].}
//' 
//' This function is a C++ implementation of \code{\link[sampling:UPMEpik2frompikw]{UPMEpik2frompikw}}.
//' More details could be find in Tille (2006).
//' @return A matrix, the joint inclusion probabilities.
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' 
//' @references 
//' Tille, Y. (2006), Sampling Algorithms, springer
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix maxentpi2(NumericVector pikr){
  
  
  double eps = 1e-6;
  int N = pikr.size();
  arma::vec pik_tmp(pikr.begin(),pikr.size(),false);
  
  // find not equl to 0 value
  arma::uvec i = arma::find(pik_tmp > eps && pik_tmp < 1-eps);
  arma::uvec i1 = arma::find(pik_tmp > 1-eps);
  
  // arma to Numeric vector 
  arma::vec pik_tmp2 = pik_tmp.elem(i);
  int N_tmp = pik_tmp2.size();
  NumericVector pik(N_tmp);
  for(int j = 0;j< N_tmp;j++){
    pik[j] = pik_tmp2[j];
  }
  
  
  NumericMatrix M2(N_tmp,N_tmp);
  NumericVector pikt = piktfrompik(pik);
  NumericVector w =  pikt/(1 - pikt);
  M2 = pik2frompik(pik,w);
  
  
  //put right value at the right place
  NumericMatrix M(N,N);
  
  
  for(int j = 0; j < i.size();j++){
    for(int ii = 0;ii < i.size();ii++){
      M(i(ii),i(j)) = M2(ii,j);
    }
    // M(_,i(j)) = M2(_,j);
  }
  for(int j = 0; j < i1.size(); j++){
    M(_,i1(j)) = pikr;
    M(i1(j),_) = pikr;
  }
  
  return(M);
}


/*** R

pik=c(0.07,0.17,0.41,0.61,0.83,0.91)

n=sum(pik)
pikt=UPMEpiktildefrompik(pik)
w=pikt/(1-pikt)
q=UPMEqfromw(w,n)
UPMEsfromq(q)

pik2frompik(pik,w)
UPMEpik2frompikw(pik,w)


pik=c(0.07,0.17,0.41,0.61,1,1,0.83,0.91)


UPmaxentropypi2(pik)
maxentpi2(pik)



pik <- inclusionprobabilities(runif(1000),100)

system.time(M1 <- round(UPmaxentropypi2(pik),4))
system.time(M2 <- round(maxentpi2(pik),4))


sum(abs(M1-M2))




*/



/*** R
library(sampling)
# pik=c(0.07,0.17,0.41,0.61,0.83,0.91)
# 
# k = 1
# set.seed(9)
pik <- sampling::inclusionprobabilities(runif(2000),500)
# w = (pik)/(1 - pik)
# qfromw(w,sum(pik))
# maxent::UPMEqfromw(w,sum(pik))
# 
# q <- qfromw(w,sum(pik))
# 
# pik2 = pik[pik != 1]
# n = sum(pik2)
# # 
# pik2
# maxent::UPMEpiktildefrompik(pik2)
# # maxent::UPMEqfromw(piktilde/(1 - piktilde),n)
# 
# pikfromq(q)
# maxent::UPMEpikfromq(q)
# 
# piktfrompik(pik2)
# maxent::UPMEpiktildefrompik(pik2)
# sum(piktfrompik(pik2))
# sum(maxent::UPMEpiktildefrompik(pik2))


system.time(s2 <- sampling::UPmaxentropy(pik))
system.time(s <- maxent(pik))

sum(s)
sum(s2)
k = k+1
.# 
# sfromq(q)
# UPMEsfromq(q)
# UPMEqfromw(w, sum(pik))




UPMEpikfromq(q)
pikfromq(q)
sum(UPMEpiktildefrompik(pik))

sum(piktfrompik(pik))




# rm(list = ls())
N=1000
n=200
pik=sampling::inclusionprobabilities(runif(N),n)
SIM=1000
PPP=rep(0,N)
for(i in 1:SIM){
  print(i)
  s <- maxent(pik)
  PPP=PPP + s
  print(sum(s))
}
PPP=PPP/SIM

t=(PPP-pik)/sqrt(pik*(1-pik)/SIM)
sum(abs(t)<1.96)/N




*/
