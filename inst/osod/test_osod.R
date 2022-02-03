
rm(list = ls())
set.seed(1)
b <- data(belgianmunicipalities)
pik <- inclusionprobabilities(belgianmunicipalities$Tot04,200)
N <- length(pik)
n <- sum(pik)
any(pik > (1-1e-7))

s <- osod(pik)
# s <- sample.speed(pik)
sum(s)




rm(list = ls())
N <- 7
pik <- c(0.5,0.5,0.3,0.1,0.6,0.7,0.3)
# pik <- pik[sample(1:N)]
s <- osod(pik)
s
SIM=10000

X <- matrix(rep(0,SIM*N),ncol = N,nrow = SIM)


for(i in 1:SIM){
  print(i)
  X[i,]=osod(pik)
  # PPP=PPP+ sample.speed(pik)
}



library("dplyr")
df <- data.frame(X)
ans = df%>%group_by_all%>%count
print("unique rows with count are")
print(ans)
View(ans)



PPP=PPP/SIM
PPP








N <- 1000
n <- 100
pik <- inclprob(runif(N),n)
s <- osod(pik)


rm(list = ls())
N <- 9000
n <- 100
pik <- inclprob(runif(N),n)

system.time(s1 <- osod(pik))
sum(s1)
system.time(s2 <- sample.speed(pik))







rm(list = ls())
pik <- c(0.4,0.5,0.8)
s <- osod(pik)
# s <- sample.speed(pik)
s
SIM=100000
PPP=rep(0,3)
for(i in 1:SIM){
  print(i)
  PPP=PPP+osod(pik)
  # PPP=PPP+ sample.speed(pik)
}
PPP=PPP/SIM
PPP


rm(list = ls())
# k = 1
# set.seed(30)
# set.seed(4)
# N=10
# n=4
# pik=inclprob(runif(N),n)
pik <- c(0.1,0.2,0.3,0.9,0.9,0.6)

# findJ2(c(0.2087,0.3130,0.9391,0.9391,0.6))
# findJ2(pik)
s <- osod(pik)
# s <- sample.speed(pik)
s
k = k+ 1
pik



rm(list = ls())
pik <- c(0.1,0.2,0.3,0.9,0.9,0.6)
SIM=1000
PPP=rep(0,length(pik))
for(i in 1:SIM){
  print(i)
  PPP=PPP+osod(pik)
  # PPP=PPP+ sample.speed(pik)
}
PPP=PPP/SIM
PPP

t=(PPP-pik)/sqrt(pik*(1-pik)/SIM)
sum(abs(t)<1.96)/5




sum(s)

rm(list = ls())
N=1000
n=200
pik=sampling::inclusionprobabilities(runif(N),n)
s <- osod(pik)
sum(s)


SIM=10000
PPP=rep(0,N)
for(i in 1:SIM){
  print(i)
  s <- osod(pik)
  PPP=PPP + s
  print(sum(s))
}
PPP=PPP/SIM

t=(PPP-pik)/sqrt(pik*(1-pik)/SIM)
sum(abs(t)<1.96)/N
