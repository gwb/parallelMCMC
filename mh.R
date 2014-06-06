
run.mh <- function(x0, niter, rproposal, dproposal, target, burnin=0){
  X <- c()
  L <- c()
  R <- c()
  X.tm1 <- x0
  for(i in seq(niter)){
    L.t <- rproposal(X.tm1)
    r <- (target(L.t) * dproposal(X.tm1, L.t)) / (target(X.tm1) * dproposal(L.t, X.tm1))

    X <- c(X, X.tm1)
    
    if(runif(1) < r){
      #X <- c(X,L.t)
      X.tm1 <- L.t
    } #else {
      #X <- c(X, X.tm1)
    #}
    
    R <- c(R,r)
    L <- c(L, L.t)
  }

  acceptance.rate <- length(unique(X))/ length(X)
  
  return(list(cbind(X[seq(burnin+1, niter)], L[seq(burnin+1,niter)], R[seq(burnin+1,niter)]), acceptance.rate))
}


run.mv.mh <- function(x0, niter, rproposal, dproposal, target, burnin=0){
  X <- NULL
  L <- NULL
  R <- c()
  X.tm1 <- x0
  acceptance.count <- 0

  for(i in seq(niter)){
    
    # proposal and acceptance
    L.t <- rproposal(X.tm1)

    if(target(X.tm1) == 0){
        warning("Density of target is null (most likely due to bad initialization), setting acceptance ratio to 1. This could have bad effect on asymptotic behavior")
        r <- 1
    } else {
      targ.L.t <- try(target(L.t), silent=T)
      if("try-error"%in% class(targ.L.t)){
        warning("target(L.t) can't be evaluated, most likely because it is far from all cluster centers. Setting acceptance ratio to 0. This could have bad effect on asymptotic behavior")
        r <- 0
      } else {
        r <- (target(L.t) * dproposal(X.tm1, L.t)) / (target(X.tm1) * dproposal(L.t, X.tm1))
      }
    }
    
    if(is.na(r)){
        save(i,x0,X,L,X.tm1, L.t, file="debug.rdata")             
    }
    
    X <- rbind(X, X.tm1)
    L <- rbind(L, L.t)

    R <- c(R, min(r,1))
    
    if(runif(1) < r) {
      acceptance.count <- acceptance.count + 1
      X.tm1 <- L.t
    }
  }
  dimnames(X) <- NULL
  dimnames(L) <- NULL

  X <- X[seq(burnin+1, niter),]
  L <- L[seq(burnin+1, niter),]
  R <- R[seq(burnin+1, niter)]
  
  return(list(X=X, L=L, R=R, acceptance.rate = acceptance.count / niter))
}


.get.constrained.target <- function(constraint.fn, original.target){
  force(constraint.fn)
  constrained.target <- function(x){
    if(constraint.fn(x)){
      return(original.target(x))
    } else {
      return(0)
    }
  }
  return(constrained.target)
}


.get.constraint.fn <- function(clust, cfn){
  force(clust)
  fn <- function(x){
    return(cfn(x) == clust)
  }
  return(fn)
}

# cfn is the cluster indicator function
# K is the number of clusters
get.constrained.targets <- function(K, cfn, original.target){
  res <- vector('list',K)
  for(i in seq(K)){
    res[[i]] <- .get.constrained.target(.get.constraint.fn(i, cfn), original.target)
    force(res)
  }
  return(res)
}
