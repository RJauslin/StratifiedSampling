context("Chauvet's method")



test_that("NO VARIABLES (only fixed size), integer in each strata",{
  N <- 100
  n <- 10
  strata <- as.matrix(rep(1:n,each = N/n))
  Xcat <- disjMatrix(strata)
  # equal
  pik <- sampling::inclusionprobastrata(strata,rep(1,n))
  X <- matrix(pik,ncol = 1)
  s <- balstrat(X,strata,pik)
  
  expect_equal(t(X/pik)%*%s,
               t(X/pik)%*%pik)
  expect_equal(sum(t(Xcat)%*%s),
               sum(t(Xcat)%*%pik))
  
  
  # unequal
  pik <- rep(sampling::inclusionprobabilities(runif(N/n),1),n)
  X <- matrix(pik,ncol = 1)
  s <- balstrat(X,strata,pik)
  
  
  expect_equal(t(X/pik)%*%s,
               t(X/pik)%*%pik)
  expect_equal(sum(t(Xcat)%*%s),
               sum(t(Xcat)%*%pik))
})


test_that("NO VARIABLES (only fixed size), real in each strata",{
  N <- 100
  n <- 10
  strata <- as.matrix(rep(1:n,each = N/n))
  Xcat <- disjMatrix(strata)
  
  # equal
  pik <- sampling::inclusionprobastrata(strata,rep(1.5,n))
  # sum(pik[strata == 1]) == 1.5
  X <- matrix(pik,ncol = 1)
  s <- balstrat(X,strata,pik)
  
  expect_equal(t(X/pik)%*%s,
               t(X/pik)%*%pik)
  expect_equal(sum(t(Xcat)%*%s),
               sum(t(Xcat)%*%pik))
  
  # unequal
  pik <- rep(sampling::inclusionprobabilities(runif(N/n),1),n)
  X <- matrix(pik,ncol = 1)
  s <- balstrat(X,strata,pik)
  
  expect_equal(t(X/pik)%*%s,
               t(X/pik)%*%pik)
  expect_equal(sum(t(Xcat)%*%s),
               sum(t(Xcat)%*%pik))
})



test_that("VARIABLES, integer in each strata",{
  N <- 100
  n <- 10
  
  p <- 4
  X <- matrix(rgamma(N*p,4,25),ncol = p)
  
  strata <- as.matrix(rep(1:n,each = N/n))
  Xcat <- disjMatrix(strata)
  
  
  # equal
  pik <- sampling::inclusionprobastrata(strata,rep(1,n))
  s <- balstrat(X,strata,pik)
  
  # sum(s)
  # t(X/pik)%*%s
  # t(X/pik)%*%pik
  
  expect_equal(abs(sum(t(Xcat)%*%s) - sum(t(Xcat)%*%pik)) < 1 + 1e-4,TRUE)
  
  
  # uneqal 
  
  pik <- rep(sampling::inclusionprobabilities(runif(N/n),1),n)
  s <- balstrat(X,strata,pik)
  
  # sum(s)
  # t(X/pik)%*%s
  # t(X/pik)%*%pik
  
  expect_equal(abs(sum(t(Xcat)%*%s) - sum(t(Xcat)%*%pik)) < 1 + 1e-4,TRUE)
  
  
})



test_that("VARIABLES, real in each strata",{
  N <- 100
  n <- 10
  
  p <- 4
  X <- matrix(rgamma(N*p,4,25),ncol = p)
  
  strata <- as.matrix(rep(1:n,each = N/n))
  Xcat <- disjMatrix(strata)
  
  
  # equal
  pik <- sampling::inclusionprobastrata(strata,rep(1.5,n))
  s <- balstrat(X,strata,pik)
  
  # sum(s)
  # t(X/pik)%*%s
  # t(X/pik)%*%pik
  
  expect_equal(abs(sum(t(Xcat)%*%s) - sum(t(Xcat)%*%pik))< 1 + 1e-4,TRUE)
  
  
  # uneqal 
  
  pik <- rep(sampling::inclusionprobabilities(runif(N/n),1.5),n)
  s <- balstrat(X,strata,pik)
  
  sum(s)
  # t(X/pik)%*%s
  # t(X/pik)%*%pik
  
  expect_equal(abs(sum(t(Xcat)%*%s) - sum(t(Xcat)%*%pik)) < 1 + 1e-4,TRUE)
  
  
})


