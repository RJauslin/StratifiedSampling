#include <RcppArmadillo.h>

using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]


void rrefArma(arma::mat& M){
  int lead = 0;
  int rowCount = M.n_rows;
  int columnCount = M.n_cols;
  double eps = 1e-11;
  int i,k;
  double temp;
  for(int r = 0; r < rowCount; r++){
    if(columnCount <= lead){
      return;
    }
    i = r;
    while(std::max(M(i,lead),-M(i,lead)) < eps ){
      M(i,lead) = 0.0;
      i = i + 1;
      if(i == rowCount){
        i = r;
        lead = lead + 1;
        if(columnCount == lead){
          return;
        }
      }
    }
    // swap rows i and r
    for(int k = 0; k < columnCount;k++){
      temp = M(i,k);
      M(i,k) = M(r,k);
      M(r,k) = temp;
    }
    // If M(r, lead) is not 0 divide row r by M(r, lead)
    if( M(r,lead) != 0.0 ){
      temp = M(r,lead);
      for(int k = 0; k < lead;k++){
        M(r,k) = 0.0;
      }
      for(int k = lead;k < columnCount;k++){
        M(r,k) = M(r,k)/temp;
      }
    }
    for(int i = 0;i < rowCount;i++){
      if( i != r ){
        // Subtract M(i, lead) multiplied by row r from row i
        temp = M(i,lead);
        for( k = 0;k < columnCount; k++){
          M(i,k) = M(i,k) - temp * M(r,k);
        }
      }
    }
    lead = lead + 1;
  }
}


/*** R
set.seed(1)
rm(list = ls())
N = 50
n = 30
p = 20
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A <- as.matrix(X/pik)
test <- t(A[1:(p+2),])
rrefArma(test)
test
*/


arma::vec osffphase(arma::vec prob, arma::mat Bm){
  int ncol = Bm.n_cols;
  int nrow = Bm.n_rows;
  
  arma::vec u(ncol,arma::fill::zeros);
  arma::uvec uset(ncol,arma::fill::zeros);
  
  double la1 = 1e+200;
  double la2 = 1e+200;
  double la, eps = 1e-9;
  int lead;
  double v, free = -1.0;
  // find nonzero vector u in Ker B (null space of B, i.e. Bu = 0)
  // with both positive and negative values
  // find reduced row echelon form of B
  rrefArma(Bm);
  
  // std::cout << Bm << std::endl;
  for(int i = (nrow-1);i >= 0; i--){
    // find lead (first nonzero entry on row) if exists
    // if no lead, i.e lead = ncol, do nothing
    // if lead, the variables after are either set or free
    // free variables are alternately set to 1 or -1
    lead = 0;
    for(int j = 0; j < ncol; j++){
      if(Bm(i,j)==0.0){
        lead++;
      }else{
        break;
      }
    }
    // lead found
    if(lead<ncol){
      v = 0.0;
      for(int j = lead+1;j < ncol;j++){
        if( uset[j] == 0 ){
          uset[j] = 1;
          free *= -1.0;
          u[j] = free;
        }
        v -= u[j]*Bm(i,j);
      }
      u[lead] = v/Bm(i,lead);
      uset[lead] = 1;
    }
  }
  // unset u[i] are free and are set to 1 or -1, can only exist at beginning
  for(int i = 0;i < ncol;i++){
    if( uset[i] == 0 ){
      free *= -1.0;
      u[i] = free;
    }else{
      break;
    }
  }
  
  // find lambda1 and lambda2
  for(int i = 0;i < ncol;i++){
    if(u[i]>0){
      la1 = std::min(la1,(1-prob[i])/u[i]);
      la2 = std::min(la2,prob[i]/u[i]);
    }
    if(u[i]<0){
      la1 = std::min(la1,-prob[i]/u[i]);
      la2 = std::min(la2,(prob[i]-1)/u[i]);
    }
  }
  // random choice of p+lambda1*u and p-lambda2*u
  if(runif(1)[0]<la2/(la1+la2)){
    la = la1;
  }else{
    la = -la2;
  }
  // update prob
  for(int i = 0;i < ncol;i++){
    prob[i] = prob[i] + la * u[i];
    if(prob[i] < eps){ prob[i] = 0; }
    if(prob[i] > 1-eps){ prob[i] = 1; }
  }
  return prob;
}

// for case where naux < p+1
void onestep_cpp_ending(arma::vec& u,arma::vec& pik,double EPS=0.0000001){
  
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

//' @title Fast flight phase of the cube method
//'
//' @description
//' 
//' This function computes the flight phase of the cube method proposed by Chauvet and Tillé (2006).
//'
//' @param Xbal A matrix of size (\eqn{N} x \eqn{p}) of auxiliary variables on which the sample must be balanced.
//' @param prob A vector of inclusion probabilities.
//' @param order if the units are reordered, Default TRUE.
//'
//'
//' @details 
//' This function implements the method proposed by (Chauvet and Tillé 2006). It recursively transforms the vector of inclusion probabilities \code{pik} into a
//' sample that respects the balancing equations. The algorithm stops when the null space of the sub-matrix \eqn{B} is empty.
//' For more information see (Chauvet and Tillé 2006).
//' 
//' 
//'
//' @return Updated vector of \code{pik} that contains 0 and 1 for unit that are rejected or selected.
//'
//' @seealso \code{\link[sampling:samplecube]{fastflightphase}}, \code{\link[BalancedSampling:cube]{cube}}.
//' 
//'
//' @examples
//' N <- 100
//' n <- 10
//' p <- 4
//' pik <- rep(n/N,N)
//' X <- cbind(pik,matrix(rgamma(N*p,4,25),ncol= p))
//' 
//' pikstar <- ffphase(X,pik) 
//' t(X/pik)%*%pikstar
//' t(X/pik)%*%pik
//' pikstar
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector ffphase(arma::mat Xbal,arma::vec prob, bool order = true){
  int N = prob.size();
  int naux = Xbal.n_cols;
  
  arma::uvec index(N);
  arma::vec p(N);
  
  int howmany;
  int k;
  for(int i = 0;i < N;i++){
    index[i]=i;
    p[i]=prob[i];
  }
  double eps = 1e-12;
  int done = 0, tempInt, howlong;
  
  // randomize order of index list
  if(order == true){
    arma::vec rnd = runif(N);
    for(int i = 0;i < N;i++){
      k = i + floor(rnd[i] * (N-i));
      tempInt = index[i];
      index[i] = index[k];
      index[k] = tempInt;
    } 
  }
  
  
  // put finished units at beginning of list
  for(int i = done;i < N;i++){
    if( p[index[i]]<eps || p[index[i]]>1-eps ){
      tempInt = index[done];
      index[done] = index[i];
      index[i] = tempInt;
      done = done + 1;
    }
  }
  
  // remaining are index from done to N-1
  while( done < N ){
    
    // find cluster of size howmany
    howmany = std::min(naux + 1,N-done);
    // if(howmany <= naux){done=N; break;}
    
    if( howmany > 1 ){
      arma::vec p_small(howmany);
      arma::vec dists(howmany); dists = 1e+20;
      arma::vec index_small(howmany);
      // arma::mat B(howmany-1,howmany);
      arma::mat B(naux,howmany);
      
      
      for(int i = 0;i < howmany; i++){
        index_small[i] = index[done+i];
        for(int j = 0;j < naux; j++){ // HERE WE CHANGED j < howmany by ----- > j < naux
          B(j,i) = Xbal(index_small[i],j)/prob[index_small[i]];
        }
        p_small[i] = p[index_small[i]];
      }
      
      if(howmany < naux + 1){
        // std::cout << "diff ending" << std::endl;
        arma::mat kern = arma::null(B);
        if(kern.empty()){
          break;
        }else{
          arma::vec u(N);
          u = kern.col(0);
          onestep_cpp_ending(u,p_small);
        }
      }else{
        p_small = osffphase(p_small,B); 
      }

      // update prob
      for(int i = 0;i < howmany;i++){
        p[index_small[i]] = p_small[i];
      }
      // update done and index
      howlong = done + howmany;
      for(int i = done;i < howlong;i++){
        if( p[index[i]]<eps || p[index[i]]>1-eps ){
          tempInt = index[done];
          index[done] = index[i];
          index[i] = tempInt;
          done = done + 1;
        }
      }
    }else{
      // max one unit left
      if(runif(1)[0] < p[index[done]]){p[index[done]]=1;}else{p[index[done]]=0;}
      done = N;
    }
  }
  
  // round
  for(int i = 0;i < N;i++){
    if( p[index[i]] > 1-eps  ){
      p[index[i]] = 1;
    }
    if( p[index[i]] < eps  ){
      p[index[i]] = 0;
    }
  }
  
  
  Rcpp::NumericVector out = Rcpp::wrap(p);
  out.attr("dim") = R_NilValue;

  return out;
}

/*** R

library(microbenchmark)
library(sampling)
rm(list = ls())
# set.seed(3)
N <-  1000
n <-  300
p <-  10
q <-  7
eps <- 1e-12
z <-  runif(N)
pik <-  inclusionprobabilities(z,n)
X <-  cbind(pik,matrix(rnorm(N*p),c(N,p)))
Z=cbind(matrix(rbinom(N*q,1,1/2),c(N,q)))
B=cbind(Z,-Z)
X <- cbind(X,B*pik)
A <- X/pik

s1 <- round(sampling::fastflightcube(X,pik),9)
s2 <- round(ffphase(X,pik),9)

*/