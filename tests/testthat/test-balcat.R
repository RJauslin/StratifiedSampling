
context("test conditional inclusion probabilities")


test_that("Simple Stratitifed design",{

  N <- 1000
  pik <- rep(0.25,N) # 8*0.25 = 2 by category
  Xcat <- as.matrix(rep(1:125,each = N/125)) # 8 by category

  Xcat_tmp <- disjMatrix(as.matrix(Xcat))
  X <- cbind(pik)
  A <- cbind(X,Xcat_tmp)/pik

  system.time(s <- balcat(X,Xcat,pik))
  expect_equal(t(A)%*%s,t(A)%*%pik)

})
