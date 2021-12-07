#' @title Sequential balanced sampling
#'
#' @description 
#' Selects at the same time a well-spread and a balanced sample using a sequential implementation.
#'
#' @param pik A vector of inclusion probabilities.
#' @param Xspread A matrix of spatial coordinates.
#' @param Xaux A matrix of auxiliary variables. The matrix must contains the \code{pik} vector to have fixed sample size.
#'
#' @details 
#' 
#' The function selects a sample using a sequential algorithm. At the same time, it respects the balancing equations (\code{Xaux}) and select a well-spread sample (\code{Xspread}). Algorithm uses a 
#' linear program to satisfy the constraints.
#'
#' @return A vector with elements equal to 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for rejected units.
#' 
#' @export
#' @importFrom stats runif
#' 
#' @examples
#' N=100
#' n <- 10
#' p=10
#' # pik=runif(N)
#' pik=inclprob(runif(N),n)
#' Xaux=array(rnorm(N*p,3,1),c(N,p))
#' 
#' Xspread <- cbind(runif(N),runif(N)) 
#' Xaux=cbind(pik,Xaux)
#' Xaux <- cbind(pik)
#' 
#' s <- balseq(pik,Xspread,Xaux)
balseq <- function(pik,Xspread,Xaux){
  
  deg = 1
  # initializing
  N <- length(pik)
  eps <- 1e-6
  
  pikInit <- pik
  index <- which(pik > eps & pik < (1-eps))
  
  n <- 0
  counter <- 1
  
  #----------- MAIN LOOP
  while(length(index) > 0){
    # cat("Step :",counter,"\n")
    # print(length(index))
    # print(n)
    
    i <- which.max(pik[index])
    
    #take distance of the considered unit 
    d <- distUnitk(Xspread,index[i],F,F)
    
    
    # modify index respect to distance
    index <- index[order(d[index])]
    
    
    l <- balseq_onestep(Xaux,pik,pikInit,index,deg)
    status <- l$status
    v = l$v
    n <- l$n
    
    unit0 <- which(pik < eps)
    unit1 <- which(pik > 1- eps)
    
    # plot(Xspread)
    # lines(Xspread[index[1:n],1],Xspread[index[1:n],2],type ="p",pch = 16,col = "cyan")
    # lines(Xspread[index[1],1],Xspread[index[1],2],type ="p",pch = 16,col = "red")
    # lines(Xspread[unit0,1],Xspread[unit0,2],type ="p",pch = 16)
    # lines(Xspread[unit1,1],Xspread[unit1,2],type ="p",pch = 16,col = "orange")
    
    # if we can no longer find solution and index is at the end of the vector then exit and return pikstar
    if(status == 1 & is.na(index[n])){
      # return(pik)
      break;
    }else{
      
      v <-  v - pmin(pik[index[2:n]],(1-pik[index[2:n]])*pik[index[1]]/(1-pik[index[1]]))
      
      if(stats::runif(1) < pik[index[1]]){
        pik[index[2:n]] <- pik[index[2:n]] - v*(1-pik[index[1]])/pik[index[1]]
        pik[index[1]] <- 1
      }else{
        pik[index[2:n]] <- pik[index[2:n]] + v
        pik[index[1]] <- 0
      }
      
      index <- which(pik > eps & pik < (1-eps))
      
    }
    counter <- counter + 1
  }
  # cat("Number of units already taked", length(which(pik > (1-eps))) ,"\n\n")
  # cat("Number of units that remains not integer",length(index),"\n\n")
  
  
  #----------- LANDING PHASE
  if(length(index) != 0){
    s <- BalancedSampling::lcube(pik,Xspread,Xaux*(pik/pikInit))
  }else{
    pik <- round(pik,10)
    s <- which(pik > (1-eps))
  }
  
  # s <- pik 
  # print(sum(pik))
  
  
  return(s)
  
}
