source("model-univariate.R")
source("performance.R") 
source("../mh.R")

mu <- 2
tau <- 0.5
X <- matrix(rtarget(1000, mu), ncol=1)
dproposal <- get.dproposal(tau)
rproposal <- get.rproposal(tau)
dtarget <- get.dtarget(mu)
dtransition <- get.dtransition(tau, mu, dproposal, dtarget)


#
get.empirical.spectral.gaps(X, 0, dtransition)
get.exact.spectral.gaps(mu, 0, dtransition, 500)

# 
get.empirical.spectral.gaps(X, 0, dproposal)
get.exact.spectral.gaps(mu, 0, dproposal, 500)



#
foo <- get.midpoint(6000, 400, 1000, rproposal, dproposal, dtarget)
Yi <- as.vector(foo[[1]])
midpoint <- foo[[2]]
clust <- foo[[4]][[1]]


min(get.empirical.spectral.gaps(Yi, midpoint, dtransition))
min(get.exact.spectral.gaps(mu, midpoint, dtransition, 300))


plot(Yi, col=clust)
abline(h=0)


# try with different midpoints
get.empirical.spectral.gaps(X, 1, dtransition)



Xv <- as.vector(X)
X1 <- matrix(Xv[Xv < 0], ncol=1)
K1 <- get.K(X1, dtransition)
K1.raw <- get.raw.K(X1, dtransition)

X2 <- matrix(Xv[Xv > 0], ncol=1)
K2 <- get.K(X2, dtransition)
K2.star <- K2[-which(diag(K2)==1), -which(diag(K2)==1)]
eigen(K2, only.values = T)$values[2]
eigen(K2.star, only.values = T)$values[2]


K1.star <- K1[-which(diag(K1)==1), -which(diag(K1)==1)]
eigen(K1, only.values = T)$values



# my version
X1 <- matrix(seq(-mu - 4*1, 0, length.out=500), ncol=1)
K1 <- get.K(X1, dproposal)
1-eigen(K1, only.values = T)$values[2]


# Aaron's
K2 <- DiscMat(-6, 0, 500, 0.5, 1)[[1]]
1-eigen(K2, only.values = T)$values[2]
