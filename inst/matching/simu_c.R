#' Title
#'
#' @param Xm 
#' @param Y 
#' @param Z 
#' @param id 
#' @param n1 
#' @param n2 
#' @param totals 
#'
#' @return
#' @export
#'
#' @examples
simu_c <- function(Xm,Y,Z,id,n1,n2,totals = FALSE){
  
  # SET UP
  y <- colnames(Y)
  z <- colnames(Z)
  N <- nrow(Xm)
  
  # samples
  pik1 <- rep(n1/N,N)
  pik2 <- rep(n2/N,N)
  
  # Xm1 <- cbind(pik1,Xm)
  # Xm2 <- cbind(pik2,Xm)
  
  # s1 <- rep(0,N)
  # s2 <- rep(0,N)
  # 
  # s1[cube(pik1,as.matrix(Xm1))] <- 1
  # s2[cube(pik2,as.matrix(Xm2))] <- 1
  # s1 <- srswor(n1,N)
  # s2 <- srswor(n2,N)
  
  s1 <- srswor(n1,N)
  comple <- which(s1 == 0)
  s2 <- rep(0,N)
  s2[comple[sample.int(length(comple),n2)]] <- 1
  
  
  X1 <- Xm[s1 == 1,]
  X2 <- Xm[s2 == 1,]
  
  
  # Y1 observed Y for S1
  # Z2 observed Z for S2
  
  Y1 <- data.frame(Y[s1 == 1,])
  colnames(Y1) <- y
  Z2 <- data.frame(Z[s2 == 1,])
  colnames(Z2) <- z
  
  id1 <- id[s1 == 1]
  id2 <- id[s2 == 1]
  
  rownames(Y1) <- id1
  rownames(Z2) <- id2
  
  d1 <- rep(N/n1,n1)
  d2 <- rep(N/n2,n2)
  
  # harmonization
  if(totals == TRUE){
    re <- harmonize(X1,d1,id1,X2,d2,id2,totals = c(N,colSums(Xm)))
  }else{
    re <- harmonize(X1,d1,id1,X2,d2,id2)
  }
  w1 = re$w1
  w2 = re$w2
  
  cat("Harmonization done \n\n")
  if(abs(sum(w1) - sum(w2)) > 1e-7){
    w2 <- w2/sum(w2)*sum(w1)
  }
  
  
  
  #----------------- RENSSEN
  
  
  fit1 <- lm(as.matrix(Y1)~X1)
  fit2 <- lm(as.matrix(Z2)~X2)
  
  
  # predict Z for S1
  Y_pred_S1 <- fit1$fitted.values
  Z_pred_S1 <- cbind(rep(1,nrow(X1)),X1)%*%fit2$coefficients
  
  # predict Y for S2
  Z_pred_S2 <- fit2$fitted.values
  Y_pred_S2 <- cbind(rep(1,nrow(X2)),X2)%*%fit1$coefficients
  
  
  m1_ren_S1 <- colSums(w1*Y_pred_S1)/sum(w1)
  m2_ren_S1 <- colSums(w1*Z_pred_S1)/sum(w1)
  
  m1_ren_S2 <- colSums(w2*Y_pred_S2)/sum(w2)
  m2_ren_S2 <- colSums(w2*Z_pred_S2)/sum(w2)
  
  
  c_ren1 <- t(w1*(Y_pred_S1 - m1_ren_S1))%*%(Z_pred_S1 - m2_ren_S1)/sum(w1)
  c_ren2 <- t(w2*(Y_pred_S2 - m1_ren_S2))%*%(Z_pred_S2 - m2_ren_S2)/sum(w2)
  
  
  #------------------- OPTIMAL
  object = otmatch(X1,id1,X2,id2,w1,w2,transport_method = "revsimplex")
  
  
  w_opt <- object[,3]
  
  y1_opt <- as.matrix(Y1[as.character(object$id1),])
  z2_opt <-  as.matrix(Z2[as.character(object$id2),])
  
  m1_opt <- colSums(w_opt*y1_opt)/sum(w_opt)
  m2_opt <- colSums(w_opt*z2_opt)/sum(w_opt)
  # m1_opt <- colMeans(y1_opt)
  # m2_opt <- colMeans(z2_opt)
  
  c_opt <- t(w_opt*(y1_opt - m1_opt))%*%(z2_opt - m2_opt)/sum(w_opt)
  
  cat("Optimal done \n\n")  
  #------------------- RANDOM
  
  out <- bsmatch(object)
  
  w_ran <- out$object[,3]
  y1_ran <- as.matrix(Y1[as.character(out$object$id1),])
  z2_ran <-  as.matrix(Z2[as.character(out$object$id2),])
  
  
  m1_ran <- colSums(w_ran*y1_ran)/sum(w_ran)
  m2_ran <- colSums(w_ran*z2_ran)/sum(w_ran)
  
  c_ran <- t(w_ran*(y1_ran - m1_ran))%*%(z2_ran - m2_ran)/sum(w_ran)
  
  
  cat("Random done \n\n")
  
  
  return(list(c_ren1 = c_ren1,
              c_ren2 = c_ren2,
              c_opt = c_opt,
              c_ran = c_ran))
  
}
