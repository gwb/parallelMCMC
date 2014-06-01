

# # # # # # # # # # # # # 
# Array helper functions
# # # # # # # # # # # # #

slice.array <- function(A, is, fn=NULL, as.str=F){
  Njs <- length(dim(A)) - length(is)
  str1 <- "A["
  str3 <- "]"
  
  if(Njs == 0){
    str2 <- paste(is, collapse=",")
  } else {
    str2.1 <- paste(is, collapse=",")
    str2.2 <- paste(rep(",", Njs), collapse="")
    str2 <- paste(c(str2.1,str2.2), collapse="")
  }

  to.eval <- paste(c(str1, str2, str3), collapse="")

  if(as.str){
    return(to.eval)
  }
  
  res <- eval(parse(text=to.eval))
  if(is.null(fn)){
    return(res)
  } else {
    return(fn(res))
  }
}


get.norm.string <- function(A,is){
  Njs <- length(dim(A)) - length(is)

  str1 <- "res["
  str3 <- "]"
  
  if(Njs == 0){
    str2 <- paste(is, collapse=",")
  } else {
    str2.1 <- paste(is, collapse=",")
    str2.2 <- paste(rep(",", Njs), collapse="")
    str2 <- paste(c(str2.1,str2.2), collapse="")
  }

  str1 <- paste(c(str1, str2, str3), collapse="")

  to.eval <- paste(str1, " <- ", str1, " / ", "norm.cst", sep="")
  
  return(to.eval)
}

# # # # # # #
# FUNCTIONS #
# # # # # # #


get.K.emp.by.index.ND <- function(discretization, X, L, alpha){
  alpha <- ifelse(alpha <= 1, alpha, 1)

  if(is.null(dim(X))){
    X <- as.matrix(X, nrow=length(X))
    L <- as.matrix(L, nrow=length(L))
  }
  
  K <- ncol(X)
  N <- nrow(X)
  
  lev <- levels(cut(X[,1], discretization))
  N.disc <- length(lev)
  rev.dic <- seq(1,N.disc)
  names(rev.dic) <- lev

  
  X.disc <- matrix(0, nrow=N, ncol=K)
  L.disc <- matrix(0, nrow=N, ncol=K)

  for(k in seq(1, K)){
    x.disc.k <- cut(X[,k], discretization)
    l.disc.k <- cut(L[,k], discretization)

    X.disc[,k] <- sapply(x.disc.k, function(a) rev.dic[a])
    L.disc[,k] <- sapply(l.disc.k, function(a) rev.dic[a])
  }

  res <- array(0, rep(N.disc, 2*K))
  norm.array <- array(0, rep(N.disc, K))

  # Careful! Normalization is not done right
  
  for(i in seq(N)){
    is <- X.disc[i,]
    js <- L.disc[i,]
    
    #norm.val <- norm.array[matrix(is, nrow=1)]
    #res.val <- res[matrix(c(is,js), nrow=1)]
    
    #res[matrix(c(is,js), nrow=1)] <- norm.val * res.val / (norm.val+1) + alpha[i]/(norm.val+1)

    res[matrix(c(is,js), nrow=1)] <-  res[matrix(c(is,js), nrow=1)] + alpha[i]
    norm.array[matrix(is, nrow=1)] <- norm.array[matrix(is, nrow=1)] + 1 #norm.val + 1
    
  }

  # This is suboptimal in terms of space. FIX IT
  coord.grid <- matrix(0, nrow=N.disc^K, ncol=K)
  for(i in seq(K)){
    coord.grid[,i] <- rep(seq(1,N.disc), each=N.disc^{i-1}, length.out=N.disc^K)
  }

  # normalization
  for(i in seq(1, N.disc^K)){
    norm.cst <- slice.array(norm.array, coord.grid[i,])
    if(norm.cst > 0){
      to.eval <- get.norm.string(res, coord.grid[i,])
      eval(parse(text=to.eval))
    }
  }

  
  for(i in seq(N.disc)){
    res[matrix(rep(i, 2*K), nrow=1)] <- 1 - slice.array(res,c(rep(i,K), rep(-i,K)), sum)
  }
  
  return(res)
}


get.K.emp.by.index.2D.tmp <- function(discretization, X, L, alpha){
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
  no.l <- colSums(res)
  for(i in seq(1, N)){
    res[i,i] <- 1 - sum(res[i,-i])
  }
  return(list(res, rbind(no.l, norm.ls) ))
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

get.K.approx.new <- function(discretization, X, L, alpha){
  alpha <- ifelse(alpha <= 1, alpha, 1)
  delta <- length(X)^2
  uniform.transition <- 1/(length(discretization) - 1)
  
  k.emp.res <- get.K.emp.by.index.new(discretization, X, L, alpha)
  k.emp <- k.emp.res[[1]]
  N <- length(discretization)-1
  k.approx <- matrix(0, nrow=N, ncol=N)
  for(i in seq(1,N)){
    for(j in seq(1,N)){
      k.approx[i,j] <- (1 - 1/delta) * k.emp[i,j] + 1/delta * 1/ (length(discretization)-1) 
    }
  }
  return(list(k.approx, k.emp.res[[2]]))
}

#attempts to fix the fact that there is equal proba to teleport everywhere
get.K.approx.new.plus <- function(discretization, X, L, alpha){
  alpha <- ifelse(alpha <= 1, alpha, 1)
  
  k.emp.res <- get.K.emp.by.index.new(discretization, X, L, alpha)
  k.emp <- k.emp.res[[1]]

  N <- length(discretization)-1
  k.approx <- matrix(0, nrow=N, ncol=N)
  for(i in seq(1,N)){
    if(k.emp[i,i] < 1){
      for(j in seq(1,N)){
        k.approx[i,j] <- k.emp[i,j]
      }
    }else{
      for(l in c(i-4, i-3, i-2, i-1, i, i+1, i+2, i+3, i+4)){
        if(l > 0 & l < N){
          k.approx[i,l] <- 1
        }
      }
        k.approx[i,] <- k.approx[i,] / sum(k.approx[i,])
    }
  }
  return(list(k.approx, k.emp.res[[2]]))
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


get.idx.fn <- function(discretization){
  fn <- function(x){
    res <- NULL

    K <- ncol(discretization) - 1
    D <- length(x)
        
    # initialization
    for(i in seq(1,D)){
      if(x[i] < discretization[i,1] || x[i] > discretization[i,K+1]){
        return(NA)
      }
    }
    
    for(i in seq(1,D)){
      for(j in seq(2, K+1)){
        if(x[i] <= discretization[i,j]){
          res <- c(res, j-1)
          break
        }
      }
    }
    return(res)
  }
  return(fn)
}

get.tsf.fn <- function(discretization){
  K <- ncol(discretization) - 1
  D <- nrow(discretization)
  fn <- function(idx){
    current.idx <- 1
    for(j in seq(1, D)){
      current.idx <- current.idx + (idx[j]-1)*K^(D-j)
    }
    return(current.idx)
  }
  return(fn)
}

# Let D be the dimension of the space, and K the number of bins
# per dimension in the discretization, and N be the number of iterations
# in the mcmc sampler.
# params: discretization is a DxK matrix
#         X is a NxD matrix
#         L is a NxD matrix
#         alpha is a vector of length N
# returns: a K^D x K^D matrix
get.K.emp.hd <- function(discretization, X, L, alpha, idx.fn, tsf.fn){
  K <- ncol(discretization) - 1
  D <- nrow(discretization)
  N <- nrow(X)
  M <- matrix(0, nrow = K^D, ncol = K^D)
  xnorm <- vector('numeric',K^D)
  
  for(i in seq(1,N)){

    idx.X <- idx.fn(X[i,])
    idx.L <- idx.fn(L[i,])

    tsf.X <- tsf.fn(idx.X)
    tsf.L <- tsf.fn(idx.L)
    
    if(is.na(tsf.X) || is.na(tsf.L)){
      warning(paste("In iteration ", i, ": X or L outside of discretization", sep=""),call.=F)
      next
    }
    
    M[tsf.X, tsf.L] <- M[tsf.X, tsf.L] + alpha[i]
    xnorm[tsf.X] <- xnorm[tsf.X] + 1
  }

  M <- M/xnorm
  M[is.na(M)] <- 0
 
  return(list(M, xnorm))
}


# note: this is tricky. In 1D, we can get away with modifying +1 and -1
# but in high D, contiguous cells might represent bins far away in space.
# FIX IT
get.K.approx.hd <- function(Kemp, xnorm, delta){
  N <- ncol(Kemp)
  for(i in seq(1, nrow(Kemp))){
    if(xnorm[i] == 0){
      if(i == 1){
        Kemp[i,i] <- 1 - (delta/2)
        Kemp[i,i+1] <- delta/2
      } else if(i == N){
        Kemp[i,i] <- 1 - (delta/2)
        Kemp[i,i-1] <- delta/2
      } else {
        Kemp[i,i] <- 1 - delta
        Kemp[i,i-1] <- delta/2
        Kemp[i,i+1] <- delta/2
      }
    }else{
      Kemp[i,] <- Kemp[i,] / sum(Kemp[i,])
    }
  }
  return(Kemp)
}

