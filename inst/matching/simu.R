#' Title
#'
#' @param Xm 
#' @param Y 
#' @param Z 
#' @param id 
#' @param n1 
#' @param n2 
#'
#' @return
#' @export
#'
#' @examples
simu <- function(Xm,Y,Z,id,n1,n2,totals = FALSE){
  
  
  # SET UP
  y <- colnames(Y)
  z <- colnames(Z)
  N <- nrow(Xm)
  YZ <- table(cbind(Y,Z))
  
  # samples
  pik1 <- rep(n1/N,N)
  pik2 <- rep(n2/N,N)

  Xm1 <- cbind(pik1,Xm)
  Xm2 <- cbind(pik2,Xm)
  # s1 <- landingRM(Xm1/pik1,ffphase(Xm1,pik1))
  # s2 <- landingRM(Xm2/pik2,ffphase(Xm2,pik2))
  
  s1 <- rep(0,N)
  s2 <- rep(0,N)
  
  s1[cube(pik1,as.matrix(Xm1))] <- 1
  s2[cube(pik2,as.matrix(Xm2))] <- 1
  
  # print(sum(s1))
  # print(sum(s2))
  
  X1 <- Xm[s1 == 1,]
  X2 <- Xm[s2 == 1,]
  
  
  Y1 <- data.frame(Y[s1 == 1,])
  colnames(Y1) <- y
  Z2 <- as.data.frame(Z[s2 == 1,])
  colnames(Z2) <- z
  
  id1 <- id[s1 == 1]
  id2 <- id[s2 == 1]
  
  rownames(Y1) <- id1
  rownames(Z2) <- id2
  
  d1 <- rep(N/n1,n1)
  d2 <- rep(N/n2,n2)
  
  # disjunctive form
  Y_dis <- sampling::disjunctive(as.matrix(Y))
  Z_dis <- sampling::disjunctive(as.matrix(Z))
  
  # Y1_dis <- sampling::disjunctive(as.matrix(Y1))
  # Z2_dis <- sampling::disjunctive(as.matrix(Z2))
  
  Y1_dis <- Y_dis[s1 ==1,]
  Z2_dis <- Z_dis[s2 ==1,]
  
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
  
  QR.1 <- qr(X1 * sqrt(w1))
  beta.yx.1 <- qr.coef(QR.1, Y1_dis * sqrt(w1))
  beta.yx.1[is.na(beta.yx.1)] <- 0
  QR.2 <- qr(X2 * sqrt(w2))
  beta.zx.2 <- qr.coef(QR.2, Z2_dis * sqrt(w2))
  beta.zx.2[is.na(beta.zx.2)] <- 0
  XX.w1 <- t(as.matrix(X1)) %*% (as.matrix(X1) * w1)
  XX.w2 <- t(as.matrix(X2)) %*% (as.matrix(X2) * w2)
  
  n12=length(intersect(id1,id2))
  gamma.p=(n1-n12)/(n1+n2-2*n12)  
  
  # gamma.p <- n1/(n1 + n2)
  XX.pool <- gamma.p * XX.w1 + (1 - gamma.p) * XX.w2
  YZ.CIA <- t(beta.yx.1) %*% XX.pool %*% beta.zx.2
  cat("Renssen done \n\n")
  
  #------------------- OPTIMAL
  object = otmatch(X1,id1,X2,id2,w1,w2)
  
  Y1_optimal <- cbind(X1[as.character(object$id1),],y = Y1[as.character(object$id1),])
  Z2_optimal <- cbind(X2[as.character(object$id2),],z = Z2[as.character(object$id2),])
  YZ_optimal <- tapply(object[,3],list(Y1_optimal$y,Z2_optimal$z),sum)
  YZ_optimal
  
  cat("Optimal done \n\n")  
  #------------------- RANDOM
  
  out <- bsmatch(object,Z2)
  
  Y1_random <- cbind(X1[as.character(out$object$id1),],y = Y1[as.character(out$object$id1),])
  Z2_random <- cbind(X2[as.character(out$object$id2),],z = Z2[as.character(out$object$id2),])
  YZ_random <- tapply(out$object$weight/out$q,list(Y1_random$y,Z2_random$z),sum)
  
  cat("Random done \n\n")
  
  return(list(YZ = YZ, YZ_ren = YZ.CIA, YZ_opt = YZ_optimal, YZ_ran = YZ_random))
  
}






























































#' Title
#'
#' @param Xm 
#' @param Y 
#' @param Z 
#' @param id 
#' @param n1 
#' @param n2 
#'
#' @return
#' @export
#'
#' @examples
simuCor <- function(Xm,Y,Z,id,n1,n2,totals = FALSE){
  
  # SET UP
  y <- colnames(Y)
  z <- colnames(Z)
  N <- nrow(Xm)
  YZ <- table(cbind(Y,Z))
  
  # samples
  s1 <- srswor(n1,N)
  s2 <- srswor(n2,N)
  
  X1 <- Xm[s1 == 1,]
  X2 <- Xm[s2 == 1,]
  
  
  Y1 <- data.frame(Y[s1 == 1,])
  colnames(Y1) <- y
  Z2 <- as.data.frame(Z[s2 == 1,])
  colnames(Z2) <- z
  
  id1 <- id[s1 == 1]
  id2 <- id[s2 == 1]
  
  rownames(Y1) <- id1
  rownames(Z2) <- id2
  
  d1 <- rep(N/n1,n1)
  d2 <- rep(N/n2,n2)
  
  # disjunctive form
  Y_dis <- sampling::disjunctive(as.matrix(Y))
  Z_dis <- sampling::disjunctive(as.matrix(Z))
  
  # Y1_dis <- sampling::disjunctive(as.matrix(Y1))
  # Z2_dis <- sampling::disjunctive(as.matrix(Z2))
  
  Y1_dis <- Y_dis[s1 ==1,]
  Z2_dis <- Z_dis[s2 ==1,]
  
  # harmonization
  # if(totals == TRUE){
  #   re <- harmonize(X1,d1,id1,X2,d2,id2,totals = c(N,colSums(Xm)))
  # }else{
  re <- harmonize(X1,d1,id1,X2,d2,id2)
  # }
  w1 = re$w1
  w2 = re$w2
  
  cat("Harmonization done \n\n")
  if(abs(sum(w1) - sum(w2)) > 1e-7){
    w2 <- w2/sum(w2)*sum(w1)
  }
  
  object = otmatch(X1,id1,X2,id2,w1,w2)
  Y1_optimal <- cbind(X1[as.character(object$id1),],y = Y1[as.character(object$id1),])
  Z2_optimal <- cbind(X2[as.character(object$id2),],z = Z2[as.character(object$id2),])
  
  mw_y <- sum(object[,3]*Y1_optimal$y)/sum(object[,3])
  mw_z <- sum(object[,3]*Z2_optimal$z)/sum(object[,3])
  
  v_y <-  t(object[,3]*(Y1_optimal$y - mw_y))%*%(Y1_optimal$y - mw_y)/sum(object[,3])
  v_z <-  t(object[,3]*(Z2_optimal$z - mw_z))%*%(Z2_optimal$z - mw_z)/sum(object[,3])
  
  cv_y_z <- (t(object[,3]*(Y1_optimal$y - mw_y))%*%((Z2_optimal$z - mw_z)))/sum(object[,3])
  
  cv_y_z/sqrt(v_y*v_z)
  cov.wt(x = cbind(Y1_optimal$y,Z2_optimal$z),wt = object[,3],cor = TRUE)$cor
  
  mean((Y$hy040n - mean(Y$hy040n))*(Z$eqIncome - mean(Z$eqIncome)))/ sqrt(mean((Y$hy040n - mean(Y$hy040n))^2)*mean((Z$eqIncome - mean(Z$eqIncome))^2))
  
}
