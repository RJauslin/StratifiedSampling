

rm(list = ls())
library(MASS)
library(BalancedSampling)
library(sampling)

Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

# test
set.seed(1)
Sigma <- Posdef(3)
N <- 10000
Xm <- mvrnorm(N,mu = c(0,0,0),Sigma)



# b1 <- matrix(c(2,-3,1),ncol = 1)
b1 <- matrix(c(0.2,-0.3,1,1.2,0.4,-0.5),ncol = 2)
Y <- Xm%*%b1 + matrix(rnorm(2*N),ncol = 2)
# Y <- Xm^2%*%b1 + matrix(rnorm(2*N),ncol = 2)


# b2 <- matrix(c(-4,1,-3),ncol = 1)
b2 <- matrix(c(-0.4,1,-0.3,-1.4,0.3,-0.6),ncol = 2)
Z <- Xm%*%b2 + matrix(rnorm(2*N),ncol = 2)
# Z <- Xm^3%*%b2 + matrix(rnorm(2*N),ncol = 2)


ALL <- cbind(Xm,Y,Z)
S <- cov(ALL)
S

# S_YX%*%(S_XX^{-1})%*%S_XZ -------------> CIA
# S_YX <- S[4,1:3]
# S_XZ <- S[1:3,5]

S_YX <- S[4:5,1:3]
S_XZ <- S[1:3,6:7]

S_YX%*%solve(S[1:3,1:3])%*%S_XZ
S
cov(Y,Z)

#------------------------------------ CIA HOLDS 

n1 <- 60
n2 <- 300
id <- seq(1,N,1)
rownames(Xm) <- id


y <- colnames(Y)
z <- colnames(Z)
N <- nrow(Xm)

s1 <- srswor(n1,N)
s2 <- srswor(n2,N)
id <- seq(1,N,1)

s1 <- srswor(n1,N)
comple <- which(s1 == 0)
s2 <- rep(0,N)
s2[comple[sample.int(length(comple),n2)]] <- 1

id1 <- id[s1 == 1]
id2 <- id[s2 == 1]
base::intersect(id1,id2)

X1 <- Xm[s1 == 1,]
X2 <- Xm[s2 == 1,]

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

######################## harmonization

re <- harmonize(X1,d1,id1,X2,d2,id2)

w1 = re$w1
w2 = re$w2

cat("Harmonization done \n\n")
if(abs(sum(w1) - sum(w2)) > 1e-7){
  w2 <- w2/sum(w2)*sum(w1)
}


######################## OTMATCH


EPS <- 1e-7
# distance
D <- proxy::dist(X1,X2,method = "Euclidean")


D_1d <- apply(D,MARGIN = 1,FUN = function(x){x})
any((D_1d[1:3000] - as.numeric(D[1,])) > 1e-7)

# to find units that are in intersection between the two sample
inter <- base::intersect(id1,id2)
inter1 <- do.call(c,lapply(inter,function(x){which(x == id1)}))
inter2 <- do.call(c,lapply(inter,function(x){which(x == id2)}))

# select the minimum weights between the two sample
wtmp <- pmin(w1[inter1],w2[inter2])

# change the weight of the intersected value to be the weight - min(w), some weights are put equal to 0.
ww1 <- w1
ww2 <- w2 
ww1[inter1] <- ww1[inter1] - wtmp
ww2[inter2] <- ww2[inter2] - wtmp

# weights that are put equal to 0, we keep only weights that are greater than 0
w01 <- ww1>EPS
w02 <- ww2>EPS

# keep track of the right indices (some weights are put to 0 and then the indices have changes)
id1_tmp <- id1[w01]
id2_tmp <- id2[w02]

# keep weights larger than 0
www1 <- ww1[w01]
www2 <- ww2[w02]

# adapt distance to optimal transport
DD <- D[w01,w02]

# optimal transport
# cat("begin transport : \n\n");start_time <- Sys.time()
# require(lpSolve)
# lp.transport(DD, "min", row.signs = rep ("==", length(www1)),
#              row.rhs = www1,
#              col.signs= rep ("==", length(www2)),
#              col.rhs = www2,
#              integers = NULL)






Amat <- function(n1,n2){
  A1 <- t(disjunctive(as.factor(rep(seq(1,n1,1),each = n2))))
  A2 <- t(disjunctive(rep(seq(1,n2,1),n1)))
  A <- rbind(A1,A2)
  return(A)
}
A <- Amat(n1,n2)
dim(A)
# image(as(A,"sparseMatrix"))


obj <- D_1d
mat <- A
rhs <- c(www1,www2)
dir = rep("==",n1+n2)


resRGLPK <- Rglpk_solve_LP(obj,
               mat,
               dir,
               rhs)

sol <- resRGLPK$solution
W <- matrix(sol,ncol = n2,byrow = T)

colSums(W)
www2
rowSums(W)
www1

res <- transport::transport.default(www1,www2,DD,method = "revsimplex")
res

W2 <- matrix(rep(0,n1*n2),ncol = n2)
for(i in 1:nrow(res)){
  W2[res$from[i],res$to[i]] <- res$mass[i]
}

colSums(W2)
www2
rowSums(W2)
www1


any((W - W2) > 1e-6)
(W-W2)[which((W - W2) > 1e-7)]







