require(kernlab)

source("models/bivariate-test.R")



get.K.hat.old <- function(X, L, nsub=500){
  Y.idx <- sample(seq(1, nrow(X)), nsub)
  Y <- X[Y.idx,]
  K.hat <- matrix(0, nrow=nsub, ncol=nsub)
  for(i in seq(1, nsub)){
    print(i)
    for(j in seq(1, nsub)){
      K.hat[i, j] <- L(Y[i,], Y[j,])
    }
  }
  return(list(K.hat=K.hat, Y=Y))
}


get.K.hat <- function(X, L, nsub=500){
  Y.idx <- sample(seq(1, nrow(X)), nsub)
  Y <- X[Y.idx,]
  get.row <- function(i) sapply(seq(nsub), function(j) L(Y[j,], Y[i,]))
  all.rows <- lapply(seq(nsub), get.row)
  return(list(K.hat=do.call('rbind', all.rows), Y=Y))
}

get.top.eigen <- function(M, num.top=NULL){
  res <- eigen(M)
  eigenvec <- res$vectors
  eigenval <- res$values
  if(!is.null(num.top)){
    eigenvec <- eigenvec[,seq(1,num.top)]
    eigenval <- eigenval[seq(1,num.top)]
  }
  return(list(vectors=eigenvec, values=eigenval))
}


subspace.proj <- function(Y, M){
  #tM <- t(M)
  #tY <- t(Y)
  #pY <- tM %*% tY
  #return(pY)
  #return(
}



mh.explorer <- get.mv.mh(3000, draw.normal.proposal, eval.normal.proposal, burnin=1000)  

init.draws <- init.exploration(dtarget, mh.explorer, matrix(c(3,3,7,-3), nrow=2,byrow=T))

X <- rtarget(1000)


K.res <- get.K.hat(X, L, nsub=200)

K.hat <- K.res$K.hat
Y <- K.res$Y

D.sq.inv <- diag(1/sqrt(rowSums(K.hat)))

K.L <- D.sq.inv %*% K.hat %*% D.sq.inv

topeigens <- get.top.eigen(K.L, 2)$vectors

normtopeigens <- t(apply(topeigens, 1, function(x) x/sqrt(sum(x^2))))

pY <- subspace.proj(t(Y), topeigens)


res <- NULL
for(j in seq(20)){
  
  nX <- rtarget(1)
  nX.dist <- sapply(seq(nrow(Y)), function(i) L(Y[i,], nX))
  norm.nX.dist <- nX.dist / sum(nX.dist)
  pnX <- c(t(norm.nX.dist) %*% normtopeigens[,1], t(norm.nX.dist) %*% normtopeigens[,2])
  norm.pnX <- pnX / sqrt(sum(pnX^2))
  res <- rbind(res, norm.pnX)
  
}


source("models/vanderwerken-1.R")
source("graphics.R")

mh.explorer <- get.mv.mh(2000, draw.normal.proposal, eval.normal.proposal, burnin=1000)  

init.points <- matrix(c(3,3,7,-3,2,7, -5,0), nrow=4, byrow=T)

init.draws <- init.exploration(dtarget, mh.explorer, init.points)

clust.res <- spectral.clustering(init.draws, 4, L, 100)

cfn <- clust.res$indicator

plot.cluster.and.target(cfn, rtarget)


foo <- get.constrained.targets(3, cfn, dtarget)
