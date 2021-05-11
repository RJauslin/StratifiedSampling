#' Title
#'
#' @param object A data.frame, output from the function \code{\link{otmatch}}.
#' @param Z2 A optional matrix, if we want to add some variables for the stratified balanced sampling step.
#'
#' @return A data.frame that contains the matching. The first two columns contain the unit identities of the two samples. The third column is the final weight. All remaining columns are the matching variables.
#' @export
#'
#' @examples
#' 
#' #--- SET UP
#' N=1000
#' p=5
#' X=array(rnorm(N*p),c(N,p))
#' EPS= 1e-9
#' 
#' n1=100
#' n2=200
#' 
#' s1=srswor(n1,N)
#' s2=srswor(n2,N)
#' 
#' 
#' id1=(1:N)[s1==1]
#' id2=(1:N)[s2==1]
#' 
#' d1=rep(N/n1,n1)
#' d2=rep(N/n2,n2)
#' 
#' X1=X[s1==1,]
#' X2=X[s2==1,]
#' 
#' #--- HARMONIZATION
#' 
#' re=harmonize(X1,d1,id1,X2,d2,id2)
#' w1=re$w1
#' w2=re$w2
#' 
#' #--- STATISTICAL MATCHING WITH OT
#' 
#' object = otmatch(X1,id1,X2,id2,w1,w2)
#' 
#' #--- 
#' 
#' out <- bsmatch(object)
#' 
#' 
#' 
#' 
#' round(colSums(object$conc$weight*object$conc),3)
#' round(colSums(cbind(w1*X1,w2*X2)),3)
#' round(colSums((out$conc$weight/out$q)*out$conc),3)
#' 
bsmatch <- function(object,
                    Z2){
  
  
  # q <- rep(0,nrow(object))
  # n1 <- nlevels(factor(object$id1))
  # l1 <- levels(factor(object$id1))
  
  q_l <- split(object$weight,f = object$id1)
  q <- as.numeric(do.call(c,lapply(q_l, function(x){x/sum(x)})))
  
  
  if(missing(Z2)){
    Z = object$weight*(object[,which(do.call(rbind,strsplit(colnames(object),"[.]"))[,1] == "X2")])  
  }else{
    Z = object$weight*(object[,which(do.call(rbind,strsplit(colnames(object),"[.]"))[,1] == "X2")])  
    Z = cbind(object$weight*disjunctive(Z2[as.character(object$id2),]),Z)
  }
  
  # strata <- cleanstrata(object$conc$id1[q < 1-EPS])
  strata <- cleanstrata(object$id1)
  # XXX <- Z[q<1-EPS,]
  XXX <- Z
  # qqq <- q[q<1-EPS]
  qqq <- q
  
  
  # strata=res[,1][q<1-EPS]
  # XXX=Z[q<1-EPS,]
  # qqq=q[q<1-EPS]
  s <- StratifiedSampling::stratifiedcube(X=XXX,strata=strata,pik=qqq)
  s <- round(s,4)
  
  # t <- tapply(s,strata,sum) 
  
  # ss <- rep(1,nrow(object$conc))
  # ss[q<1-EPS] <- s
  ss <- s
  
  out <- object
  # out$res <- out$res[ss == 1,]
  # out <- out[ss == 1,]
  # out$q <- q[ss == 1]
  out <- list(object = object[ss == 1,],q = q[ss == 1])
  
  
  return(out)
  
}




