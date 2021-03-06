% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/harmonize.R
\name{harmonize}
\alias{harmonize}
\title{Harmonization by calibration}
\usage{
harmonize(X1, d1, id1, X2, d2, id2, totals)
}
\arguments{
\item{X1}{A matrix, the matching variables of sample 1.}

\item{d1}{A numeric vector that contains the initial weights of the sample 1.}

\item{id1}{A character or numeric vector that contains the labels of the units in sample 1.}

\item{X2}{A matrix, the matching variables of sample 2.}

\item{d2}{A numeric vector that contains the initial weights of the sample 1.}

\item{id2}{A character or numeric vector that contains the labels of the units in sample 2.}

\item{totals}{An optional numeric vector that contains the totals of the matching variables.}
}
\value{
A list of two vectors, the new weights of sample 1 (respectively new weights of sample 2).
}
\description{
This function harmonize the two weight schemes such that the totals are equal.
}
\details{
All details of the method can be seen in the manuscript: Raphaël Jauslin and Yves Tillé (2021) <arXiv:>.
}
\examples{

#--- SET UP

N = 1000
p = 5
X = array(rnorm(N*p),c(N,p))

n1=100
n2=200

s1 = sampling::srswor(n1,N)
s2 = sampling::srswor(n2,N)

id1=(1:N)[s1==1]
id2=(1:N)[s2==1]

d1=rep(N/n1,n1)
d2=rep(N/n2,n2)

X1 = X[s1==1,]
X2 = X[s2==1,]

re <- harmonize(X1,d1,id1,X2,d2,id2)

colSums(re$w1*X1)
colSums(re$w2*X2)

#---  if the true totals is known

totals <- c(N,colSums(X))
re <- harmonize(X1,d1,id1,X2,d2,id2,totals)

colSums(re$w1*X1)
colSums(re$w2*X2)
colSums(X)

}
