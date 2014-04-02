

# # # # # # #
# FUNCTIONS #
# # # # # # #


get.K.emp.by.index.2D <- function(discretization, X, L, alpha){
  alpha <- ifelse(alpha <= 1, alpha, 1)

  lev <- cut(X[,1], discretization)

  X.disc <- apply(X, 2, function(mycol) cut(mycol, discretization))
  L.disc <- apply(L, 2, function(mycol) cut(mycol, discretization))

  N <- length(lev)
  
  rev.disc <- array(0, c(N, N, 2))
  for(i in seq(1, N)){
    for(j in seq(1, N)){
      rev.disc[i,j,] <- c(i,j)
    }
  }
  dimnames(rev.disc) <- list(lev, lev)

  res <- array(0, c(N, N, N, N))
  norm.ls <- matrix(0, nrow=N, ncol=N)

  for(k in seq(length(X))){
    is <- rev.disc[X.disc[k,1], X.disc[k,2]]
    js <- rev.disc[L.disc[k,1], L.disc[k,2]]
    al <- alpha[k]
    res[is[1],is[2], js[1], js[2]] <- res[is[1],is[2], js[1], js[2]] + al
    norm.ls[is[1], is[2]] <- norm.ls[is[1], is[2]] + 1
  }
  
  
}

get.K.emp.by.index.new <- function(discretization, X, L, alpha){
  alpha <- ifelse(alpha <= 1, alpha, 1)
  X.disc <- cut(X, discretization)
  L.disc <- cut(L, discretization)

  rev.disc <- seq(length(X.disc))
  names(rev.disc) <- levels(X.disc)

  N <- length(levels(X.disc))
  
  res <- matrix(0, nrow=N, ncol=N)
  norm.ls <- vector('numeric', length=N)

  for(k in seq(length(X))){
    i <- rev.disc[X.disc[k]]
    j <- rev.disc[L.disc[k]]
    al <- alpha[k]
    res[i,j] <- res[i,j] + al
    norm.ls[i] <- norm.ls[i] + 1
  }

  res <- res / norm.ls
  res <- ifelse(is.na(res), 0, res)
  for(i in seq(1, N)){
    res[i,i] <- 1 - sum(res[i,-i])
  }
  return(res)
}


get.K.emp.by.index <- function(discretization, X, L, alpha){
  alpha <- ifelse(alpha <= 1, alpha, 1)
  X.disc <- cut(X, discretization)
  L.disc <- cut(L, discretization)
  k.emp.tmp <- function(i,j){
    # i and j represent indices in the partition
    # alpha is the acceptance function (aka acceptance ratio)

    i.int <- levels(X.disc)[i]
    j.int <- levels(X.disc)[j]
    
    norm.cst <- sum(X.disc == i.int)
    transit <- sum( (X.disc == i.int & L.disc == j.int) * alpha )
    
    if(norm.cst == 0){
      return(0)
    } else {
      return( transit / norm.cst )
    }
  }

  N <- length(discretization) - 1
  
  kmat <- matrix(0, nrow=N, ncol=N)
  for(i in seq(1,N)){
    for(j in seq(1,N)){
      kmat[i,j] <- k.emp.tmp(i,j)
    }
  }

  for(i in seq(1,N)){
    kmat[i,i] <- 1 - sum(kmat[i,-i])
  }
  return(kmat)
}



get.K.approx <- function(discretization, X, L, alpha){
  alpha <- ifelse(alpha <= 1, alpha, 1)
  delta <- length(X)^2
  uniform.transition <- 1/(length(discretization) - 1)
  
  k.emp <- get.K.emp.by.index(discretization, X, L, alpha)
  N <- length(discretization)-1
  k.approx <- matrix(0, nrow=N, ncol=N)
  for(i in seq(1,N)){
    for(j in seq(1,N)){
      k.approx[i,j] <- (1 - 1/delta) * k.emp[i,j] + 1/delta * 1/ (length(discretization)-1) 
    }
  }
  return(k.approx)
}

get.K.approx <- function(discretization, X, L, alpha){
  alpha <- ifelse(alpha <= 1, alpha, 1)
  delta <- length(X)^2
  uniform.transition <- 1/(length(discretization) - 1)
  
  k.emp <- get.K.emp.by.index.new(discretization, X, L, alpha)
  N <- length(discretization)-1
  k.approx <- matrix(0, nrow=N, ncol=N)
  for(i in seq(1,N)){
    for(j in seq(1,N)){
      k.approx[i,j] <- (1 - 1/delta) * k.emp[i,j] + 1/delta * 1/ (length(discretization)-1) 
    }
  }
  return(k.approx)
}


# used to have a sense for wether the empirical transition makes any sense
simple.riemann <- function(breaks, fn){
  step <- diff(breaks[1:2])
  return(sum(sapply(breaks, fn) * step))
}


# checking wether the empirical transition makes any sense
# tst.fn <- function(x){
#   return(pnorm(1.06, x, 0.1) - pnorm(1.01, x, 0.1))
# }


