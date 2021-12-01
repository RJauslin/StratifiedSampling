
rm(list = ls())

## initialization

N=100
n <- 10
p=10
# pik=runif(N)
pik=inclusionprobabilities(runif(N),n)
Xaux=array(rnorm(N*p,3,1),c(N,p))

Xspread <- cbind(runif(N),runif(N))
Xaux=cbind(pik,Xaux)
# Xaux <- cbind(pik)


# one step 




N <- length(pik)
eps <- 1e-6






J <- 4
pik1 <- pik[1:J]
pik2 <- pik[2:J]


# X <- matrix(Xaux[1:J,],ncol = ncol(Xaux))
# X <- X*(pik1/pikInit[index[1:n]])


x1 <- Xaux[1,]
ck <- c()
for(k in 2:J){
  ck <- c(ck,(pik[1]/(pik[k]*(1- 2*pi[1])) )*t(1/x1)%*%Xaux[k,])
}



# bmin <- -pik[1:J]
# bmax <- 1 - pik[1:J]
bmin <- -pmin(pik1,(1 - pik1) * pik1[1]/(1-pik1[1]))
bmax <- pmin(1 - pik1,pik1 * pik1[1]/(1-pik1[1]))

mat <- rbind(ck,diag(rep(1,J-1)))
rhs <- c(bmax-bmin)
dir = rep("<=",J)



out <- Rglpk::Rglpk_solve_LP(obj = ck,
                      mat = mat,
                      dir = dir,
                      rhs = c(bmax-bmin),
                      max = TRUE
)

pikInit <- pik

# option 1
pik[1] <- pik[1] + (out$optimum - pik[1])
pik[2:J] <- pik[2:J] + (out$solution - pik[2:J])

pikInit[1:J]
pik[1:J]

# option 2



pikInit <- pik
index <- which(pik > eps & pik < (1-eps))


counter <- 1

i <- which.max(pik[index])
d <- distUnitk(Xspread,index[i],F,F)
index <- index[order(d[index])]



# n <- ncol(Xaux)- 1
# p <- ncol(Xaux) -1

J <- 15
p <- ncol(Xaux)

status <- 1
sol <- 0
# loop on the linear programming
  
  
pik1 <- pik[index[1:n]]
pik2 <- pik[index[2:n]]
  
bmin <- -pmin(pik2,(1 - pik2) * pik1[1]/(1-pik1[1]))
bmax <- pmin(1 - pik2,pik2 * pik1[1]/(1-pik1[1]))
  
# X = matrix(Xaux[index[1:n],],ncol = ncol(Xaux))
# X <- Xaux[index[1:n],]
X <- matrix(Xaux[index[1:n],],ncol = ncol(Xaux))
X <- X*(pik1/pikInit[index[1:n]])



V <- Rglpk::Rglpk_solve_LP(obj = ((n-1):1),
                           mat = rbind(t(X[-1,]/pik2),diag(rep(1,n-1))),
                           dir = c(rep("==",p+1),rep("<=",n-1)),
                           rhs = c(X[1,] + c(t(X[-1,]/pik2)%*%pmin(pik2,(1-pik2)*pik1[1]/(1-pik1[1]))),
                                 bmax-bmin),
                           max = TRUE
)
sol <- V$solution
  
status <- V$status






# 
# 
# l <- balseq_onestep(Xaux,pik,pikInit,index,deg)
# status <- l$status
# v = l$v
# n <- l$n
#   
# unit0 <- which(pik < eps)
# unit1 <- which(pik > 1- eps)
# 


