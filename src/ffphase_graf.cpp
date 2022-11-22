#include <RcppArmadillo.h>

using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
//' @title reduced row echelon form arma implementation
//'
//'
//' @param M matrix 
//'
//' @return NULL (transform matrix)
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
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


// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//' @description
//' description
//'
//'
//' @param prob inclusion probabilities
//' @param Bm matrix of auxiliary variables
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


// [[Rcpp::depends(RcppArmadillo)]]
//' @title Fast Flight phase
//'
//'
//' @description
//' 
//' Modified version of \code{\link[BalancedSampling:flightphase]{flightphase}}
//'
//' @param prob vector of inclusion probabilities of size N.
//' @param Xbal Matrix of auxiliary variables of dimension N x p
//' @param order if reordering at first step, Default TRUE.
//' @param redux if the matrix should be reduced. Default FALSE.
//'
//' @details
//'
//' details
//'
//' @return a sample with at most p value not update to 0 or 1. 
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::vec ffphase_graf(arma::vec prob, arma::mat Xbal, bool order = true, bool redux = false){
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
    if(howmany <= naux){done=N; break;}
    
    if( howmany > 1 ){
      arma::vec p_small(howmany);
      arma::vec dists(howmany); dists = 1e+20;
      arma::vec index_small(howmany);
      arma::mat B(howmany-1,howmany);
      // arma::mat B(naux,howmany);
      
      
      for(int i = 0;i < howmany; i++){
        index_small[i] = index[done+i];
        for(int j = 0;j < howmany-1;j++){ // HERE WE CHANGED j < howmany by ----- > j < naux
          B(j,i) = Xbal(index_small[i],j)/prob[index_small[i]];
        }
        p_small[i] = p[index_small[i]];
      }
      p_small = osffphase(p_small,B); 

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
  return p;
}

/*** R


library(sampling)
rm(list = ls())
N = 3000
n = 200
p = 10
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))

ffphase_graf(pik,X)

library("microbenchmark")

micro_timings = microbenchmark(ffphase(X,pik),
                               ffphase_cpp(X,pik),
                               ffphase_graf(pik,X),
                               BalancedSampling::flightphase(pik,X),
                               times = 30)


plot(micro_timings)
micro_timings











rm(list = ls())
# set.seed(3)
N <-  1000
n <-  300
p <-  1
q <-  7
eps <- 1e-12
z <-  runif(N)
pik <-  inclusionprobabilities(z,n)
X <-  cbind(pik,matrix(rnorm(N*p),c(N,p)))
Z=cbind(matrix(rbinom(N*q,1,1/2),c(N,q)))
B=cbind(Z,-Z)
X <- cbind(X,B*pik)
A <- X/pik
pikfastflightcube <- round(fastflightcube(X,pik),9)
pikffphase <- round(SamplingC::ffphase(pik,X),9)
pikBal <-  round(SamplingC::flightphase(pik,X),9)
dim(X)
length(which(pikfastflightcube > eps & pikfastflightcube < (1-eps)))
length(which(pikBal > eps & pikBal < (1-eps)))
length(which(pikffphase > eps & pikffphase < (1-eps)))
t(A)%*%pikfastflightcube
t(A)%*%pikBal
t(A)%*%pikffphase
##########################################################################
rm(list = ls())
eps = 1e-12
N = 5000
n = 800
p = 40
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
system.time(test1 <- SamplingC::flightphase(pik,X))
system.time(test2 <- ffphase(pik,X))
length(which(test1 > eps & test1 < (1-eps)))
length(which(test2 > eps & test2 < (1-eps)))
##############################################################################
rm(list = ls())
set.seed(1)
eps <- 1e-13
library(Matrix)
N <- 200
n1 <- floor(N/3)
n2 <- floor(N/5)
n3 <- floor(N/7)
Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),n1),
                sampling::inclusionprobabilities(runif(N),n2),
                sampling::inclusionprobabilities(runif(N),n3)),ncol = 3)
X <- PM(Pik)$PM
X <- cbind(rep(1,nrow(X)),X)
pik <- PM(Pik)$P
image(as(X,"sparseMatrix"))
A <- X/pik
system.time(test <- SamplingC::ffphase(pik,X,order = TRUE,redux = FALSE)) # order changed and no reduc.
# t(A)%*%test
system.time(test <- SamplingC::ffphase(pik,X,order = TRUE,redux = TRUE)) # order changed with reduc.
# t(A)%*%test
system.time(test <- SamplingC::ffphase(pik,X,order = FALSE,redux = TRUE)) # no change in the order and with reduc.
# t(A)%*%test
# t(A)%*%pik
system.time(test <- sampling::fastflightcube(X,pik,order = 2))
system.time(test <- BalancedSampling::flightphase(pik,X))
# system.time(test <- flightphase_arma(X,pik))
# system.time(test <- flightphase_arma2(X,pik))
t(A)%*%pik
t(A)%*%test
t(A)%*%pik # correct
t(A)%*%pikCube01
length(which(test > eps & test < (1-eps)))
length(which(test > eps & test < (1-eps)))
test <- matrix(c(1:6,rep(0,6*4)),ncol = 6,nrow = 5,byrow = T)
t <- matrix(rnorm(5*6),ncol = 6,nrow = 5)
rrefBal(t)
ukern(t)
library(MASS)
Null(t(test))
rm(list = ls())
N = 50
n = 30
p = 2
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
flightphaseSPOT(pik,X)
*/