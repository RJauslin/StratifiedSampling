#' Title
#'
#' @param X 
#' @param strata 
#' @param pik 
#' @param comment 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' rm(list = ls())
#' N <- 10000
#' n <- c(25,50,100,250)
#' p <- 3
#' X_tmp <- matrix(rgamma(N*p,4,25),ncol = p)
#' 
#' l <- list()
#' sizeStrata <- c(1,1.5)
#' step <- 1
#' 
#' i = 4
#' j = 1
#' strata <- as.matrix(rep(1:n[i],each = N/n[i]))
#' pik <- inclusionprobastrata(strata,rep(sizeStrata[j],n[i]))
#' X <- cbind(pik,X_tmp)
#' comment = TRUE
#' 
#' balancedstratification2(X,strata,pik,comment=TRUE))[3]
#' 
balancedstratification2 <- function (X, strata, pik) 
{
  
  
  
  ##----------------------------------------------------------------
  ##                          Initializing                         -
  ##----------------------------------------------------------------
  
  
  EPS = 1e-8
  strata = sampling::cleanstrata(strata)
  H = max(strata)
  N = dim(X)[1]
  pikstar = rep(0, times = N)
  
  
  ##----------------------------------------------------------------
  ##                            Step 1                             -
  ##----------------------------------------------------------------
  
  
  for (h in 1:H) {
    pikstar[strata == h] = ffphase(cbind(pik[strata == h],X[strata == h, ]), pik[strata == h])
  }
  
  
  ##----------------------------------------------------------------
  ##                            Step 2                             -
  ##----------------------------------------------------------------
  
  
  XN = cbind(sampling::disjunctive(strata) * pik, X)/pik * pikstar
  if (is.null(colnames(X)) == FALSE) 
    colnames(XN) <- c(paste("Stratum", 1:H, sep = ""), 
                      colnames(X))
  
  i <- which(pikstar > EPS & pikstar < (1-EPS))
  pikstar[i] = ffphase(as.matrix(cbind(X[i,],XN[i,])), pikstar[i])
  
  
  
  ##----------------------------------------------------------------
  ##                            Step 3                             -
  ##----------------------------------------------------------------
  
  
  i <- which(pikstar > EPS & pikstar < (1-EPS))
  pikstar[i] <- landingRM(as.matrix(cbind(X[i,],XN[i,])),
                          pikstar[i],
                          pik[i])
  pikstar <- round(pikstar,8)
  
  
  if (comment) {
    A_tmp <- as.matrix(cbind(sampling::disjunctive(strata)* pik,X))
    A = A_tmp[pik > EPS, ]/pik[pik > EPS]
    TOT = t(A) %*% pik[pik > EPS]
    EST = t(A) %*% pikstar[pik > EPS]
    DEV = 100 * (EST - TOT)/TOT
    cat("\n\nQUALITY OF BALANCING\n")
    if (is.null(colnames(X))){
      Vn = as.character(1:length(TOT))
    }else{
      Vn = colnames(A)
    }
    for (i in 1:length(TOT)) if (Vn[i] == "") 
      Vn[i] = as.character(i)
    d = data.frame(TOTALS = c(TOT), HorvitzThompson_estimators = c(EST), 
                   Relative_deviation = c(DEV))
    rownames(d) <- Vn
    print(d)
  }
  return(pikstar)
}