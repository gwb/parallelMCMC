require(mvtnorm)
require(clusterGeneration)


C <- 4
d <- 5

sig.ls <- vector('list', C)
u.ls <- vector('list', C)

w.ls <- c(0.1, 0.2, 0.3, 0.4)

sig.prop <- 0.2*diag(d)

for(i in seq(C)){
  sig.ls[[i]] <- genPositiveDefMat(d, covMethod="onion", rangeVar=c(0.2,1))$Sigma
  u.ls[[i]] <- runif(d, min=-10, max=10)
}

dtarget <- function(x){
  res <- 0
  for(i in seq(C)){
    res <- res + w.ls[i] * dmvnorm(x, u.ls[[i]], sig.ls[[i]])
  }
  return(res)
}

rtarget <- function(n){
  res <- NULL
  for(j in seq(n)){
    i <- sample(seq(C), 1, prob=w.ls)
    res <- rbind(res, rmvnorm(1, u.ls[[i]], sig.ls[[i]]))
  }
  return(res)
}

draw.normal.proposal <- function(x){
  return(as.vector(rmvnorm(1, x, sig.prop)))
}

eval.normal.proposal <- function(new.state, old.state){
  return(dmvnorm(new.state, mean=old.state, sigma=sig.prop))
}

# make sure the function never returns NA (take care of over/underflows)
L <- function(new.state, old.state){
  r.num <- dtarget(new.state) * eval.normal.proposal(old.state, new.state)
  r.denom <- dtarget(old.state) * eval.normal.proposal(new.state, old.state)
  a <- min(1, r.num/r.denom)
  res <- eval.normal.proposal(new.state, old.state) * a
  if(is.na(res)){
    return(0)
  } else {
    return(res)
  }
}
