require(mvtnorm)

# # # # # # # # # # # # # # # # # # # # # # 
# Model (same target as vanderwerken 2013)
# # # # # # # # # # # # # # # # # # # # # # 

sig1 <- matrix(c(1,0.2,0.2,1), byrow=T, nrow=2)
sig2 <- matrix(c(2,-0.5, -0.5, 0.5), byrow=T,nrow=2)
sig3 <- matrix(c(1.3, 0.3, 0.3, 0.4), byrow=T, nrow=2)
sig4 <- matrix(c(1,1,1,2.5), byrow=T, nrow=2)

u1 <- c(3,3)
u2 <- c(7,-3)
u3 <- c(2,7)
u4 <- c(-5,0)

sig.ls <- list(sig1, sig2, sig3, sig4)
u.ls <- list(u1, u2, u3, u4)
w.ls <- c(0.02, 0.2, 0.2, 0.58)


sig.prop <- matrix(c(0.2,0,0,0.2), byrow=T, nrow=2)


dtarget <- function(x){
  res <- 0
  for(i in seq(4)){
    res <- res + w.ls[i] * dmvnorm(x, u.ls[[i]], sig.ls[[i]])
  }
  return(res)
}

rtarget <- function(n){
  res <- NULL
  for(j in seq(n)){
    i <- sample(c(1,2,3,4), 1, prob=w.ls)
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
