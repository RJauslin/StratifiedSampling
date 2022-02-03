rm(list = ls())
N <- 10000
n <- 1000
# x1 <- rgamma(N,4,25)
# x2 <- rgamma(N,4,25)

strata <- as.matrix(rep(1:n,each = N/n))
Xcat <- disjMatrix(strata)
pik <- inclusionprobastrata(strata,rep(1.5,n))
# X <- as.matrix(cbind(matrix(c(x1,x2),ncol = 2)))
p <- 10
X <- as.matrix(matrix(rgamma(N*p,4,25),ncol = p))



system.time(s <- stratifiedcube(X,strata,pik))
system.time(s3 <- fbs(X,strata,pik))

t(X/pik)%*%s
t(X/pik)%*%pik
t(X/pik)%*%s3

t(Xcat_tmp)%*%s
t(Xcat_tmp)%*%s3

sum(s3)




########################## NOT EQUAL TO INTEGER INSIDE STRATA


rm(list = ls())
N <- 10000
n <- 1000
# x1 <- rgamma(N,4,25)
# x2 <- rgamma(N,4,25)

Xcat <- as.matrix(rep(1:n,each = N/n))
Xcat_tmp <- disjMatrix(Xcat)


# pik <- inclusionprobastrata(Xcat,rep(1,n))
pik <- rep(0.25,N)
sum(pik[Xcat == 1])


# X <- as.matrix(cbind(matrix(c(x1,x2),ncol = 2)))
p <- 10
X <- as.matrix(matrix(rgamma(N*p,4,25),ncol = p))

# X <- as.matrix(pik)

system.time(s <- balcat(X,Xcat,pik))
system.time(s3 <- fbstratification(X,Xcat,pik))

t(X/pik)%*%s
t(X/pik)%*%pik
t(X/pik)%*%s3

t(Xcat_tmp)%*%s
t(Xcat_tmp)%*%pik

t(Xcat_tmp)%*%s3

sum(s3)


