simu2 <- function(Xm,Y,Z,id,n1,n2,totals = FALSE){
  
  # SET UP
  y <- colnames(Y)
  z <- colnames(Z)
  N <- nrow(Xm)
  # res_true <- tapply(YZ$income,list(YZ$ecostat),mean)
 
  
  
  # balanced sampling 
  pik1 <- rep(n1/N,N)
  pik2 <- rep(n2/N,N)
  Xm1 <- cbind(pik1,Xm)
  Xm2 <- cbind(pik2,Xm)
  s1 <- rep(0,N)
  s2 <- rep(0,N)
  s1[cube(pik1,as.matrix(Xm1))] <- 1
  s2[cube(pik2,as.matrix(Xm2))] <- 1
  
  # srsrwor
  
  # s1 <- srswor(n1,N)
  # s2 <- srswor(n2,N)
  
  # selected samples
  
  
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
  Y1_dis <- Y_dis[s1 ==1,]
  
  
  re <- harmonize(X1,d1,id1,X2,d2,id2,totals = c(N,colSums(Xm)))
  
  w1 = re$w1
  w2 = re$w2
  
  # cat("Harmonization done \n\n")
  if(abs(sum(w1) - sum(w2)) > 1e-7){
    w2 <- w2/sum(w2)*sum(w1)
  }
  
  #----------------- optimal transport
  
  object <- otmatch(X1,id1,X2,id2,w1,w2)
  count_opt <- tapply(object$weight*Z2[as.character(object$id2),],list(Y1[as.character(object$id1),]),sum)
  m_opt <- tapply(object$weight,list(Y1[as.character(object$id1),]),sum)
  
  res_opt <- count_opt/m_opt
  
  #------------------- balanced sampling match
  
  out <- bsmatch(object)
  
  count_ran <- tapply( (out$object$weight/out$q)*Z2[as.character(out$object$id2),],list(Y1[as.character(out$object$id1),]),sum)
  m_ran <- tapply(out$object$weight/out$q,list(Y1[as.character(out$object$id1),]),sum)
  res_ran <- count_ran/m_ran
  
  #----------------- RENSSEN
  
  QR.1 <- qr(X1 * sqrt(w1))
  beta.yx.1 <- qr.coef(QR.1, Y1_dis * sqrt(w1))
  beta.yx.1[is.na(beta.yx.1)] <- 0
  QR.2 <- qr(X2 * sqrt(w2))
  beta.zx.2 <- qr.coef(QR.2, Z2 * sqrt(w2))
  beta.zx.2[is.na(beta.zx.2)] <- 0
  XX.w1 <- t(as.matrix(X1)) %*% (as.matrix(X1) * w1)
  XX.w2 <- t(as.matrix(X2)) %*% (as.matrix(X2) * w2)
  
  n12=length(intersect(id1,id2))
  gamma.p=(n1-n12)/(n1+n2-2*n12)  
  
  # gamma.p <- n1/(n1 + n2)
  XX.pool <- gamma.p * XX.w1 + (1 - gamma.p) * XX.w2
  YZ.CIA <- t(beta.yx.1) %*% XX.pool %*% beta.zx.2
  res_ren <- as.vector(YZ.CIA)/tapply(w1,list(Y1$ecostat),sum)  
  
  
  #----------------- return
  return(list(res_opt = res_opt,
              res_ran = res_ran,
              res_ren = res_ren))
  
}



simu_strata <- function(Xm,Y,Z,id,pik1,pik2,strata,totals = FALSE){
  
  # SET UP
  y <- colnames(Y)
  z <- colnames(Z)
  N <- nrow(Xm)
  # res_true <- tapply(YZ$income,list(YZ$ecostat),mean)
  
  n1 <- sum(pik1)
  n2 <- sum(pik2)
  
  # stratified balanced sampling
  
  
  
  s1 <- stratifiedcube(X = as.matrix(Xm),
                      strata = strata,
                      pik = pik1)
  s2 <- stratifiedcube(X = as.matrix(Xm),
                       strata = strata,
                       pik = pik2)
  
  
  
  # selected samples
  
  
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
  
  
  d1 <- 1/pik1[s1 == 1]
  d2 <- 1/pik2[s2 == 1]
  
  # d1 <- rep(N/n1,n1)
  # d2 <- rep(N/n2,n2)
  
  # disjunctive form
  Y_dis <- sampling::disjunctive(as.matrix(Y))
  Y1_dis <- Y_dis[s1 ==1,]
  
  
  re <- harmonize(X1,d1,id1,X2,d2,id2,totals = c(N,colSums(Xm)))
  
  w1 = re$w1
  w2 = re$w2
  
  # cat("Harmonization done \n\n")
  if(abs(sum(w1) - sum(w2)) > 1e-7){
    w2 <- w2/sum(w2)*sum(w1)
  }
  
  #----------------- optimal transport
  
  object <- otmatch(X1,id1,X2,id2,w1,w2)
  count_opt <- tapply(object$weight*Z2[as.character(object$id2),],list(Y1[as.character(object$id1),]),sum)
  m_opt <- tapply(object$weight,list(Y1[as.character(object$id1),]),sum)
  
  res_opt <- count_opt/m_opt
  
  #------------------- balanced sampling match
  
  out <- bsmatch(object)
  
  count_ran <- tapply( (out$object$weight/out$q)*Z2[as.character(out$object$id2),],list(Y1[as.character(out$object$id1),]),sum)
  m_ran <- tapply(out$object$weight/out$q,list(Y1[as.character(out$object$id1),]),sum)
  res_ran <- count_ran/m_ran
  
  #----------------- RENSSEN
  
  QR.1 <- qr(X1 * sqrt(w1))
  beta.yx.1 <- qr.coef(QR.1, Y1_dis * sqrt(w1))
  beta.yx.1[is.na(beta.yx.1)] <- 0
  QR.2 <- qr(X2 * sqrt(w2))
  beta.zx.2 <- qr.coef(QR.2, Z2 * sqrt(w2))
  beta.zx.2[is.na(beta.zx.2)] <- 0
  XX.w1 <- t(as.matrix(X1)) %*% (as.matrix(X1) * w1)
  XX.w2 <- t(as.matrix(X2)) %*% (as.matrix(X2) * w2)
  
  n12=length(intersect(id1,id2))
  gamma.p=(n1-n12)/(n1+n2-2*n12)  
  
  # gamma.p <- n1/(n1 + n2)
  XX.pool <- gamma.p * XX.w1 + (1 - gamma.p) * XX.w2
  YZ.CIA <- t(beta.yx.1) %*% XX.pool %*% beta.zx.2
  res_ren <- as.vector(YZ.CIA)/tapply(w1,list(Y1$ecostat),sum)  
  
  
  #----------------- return
  return(list(res_opt = res_opt,
              res_ran = res_ran,
              res_ren = res_ren))
  
}

