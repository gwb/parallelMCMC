require(mvtnorm)
require(clusterGeneration)


C <- 4
d <- 5

sig.ls <- vector('list', C)
u.ls <- vector('list', C)

w.ls <- c(0.1, 0.2, 0.3, 0.4)

sig.prop <- 0.2*diag(d)

for(i in seq(C)){
  sig.ls[[i]] <- diag(d)
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

rtarget.clust <- function(n){
  res <- NULL
  clust <- NULL
  for(j in seq(n)){
    i <- sample(seq(C), 1, prob=w.ls)
    res <- rbind(res, rmvnorm(1, u.ls[[i]], sig.ls[[i]]))
    clust <- c(clust, i)
  }
  return(list(res, clust))
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



# this implementations leverages the fact that the
# proposal is symmetric. I'm reluctant because it means that
# the algorithm is still slow in general, but for prototyping
# purposes, I need all the tricks I can get
L.vec.sym <- function(new.states, old.state){
  eval.proposals <- eval.normal.proposal(new.states,old.state)
  r.num <- dtarget(new.states) * eval.proposals
  r.denom <- dtarget(old.state) * eval.proposals
  a <- pmin(1, r.num/r.denom)
  res <- eval.proposals * a
  return(ifelse(is.na(res), 0, res))
}
