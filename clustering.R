

cluster.kmeans <- function(X,K,iter.max=50,nstart=20){
  # each row of X should be a point. so the number of columns of X should be the dimension of the space
  res <- kmeans(X,K,iter.max=iter.max,nstart=nstart)
  return(res) #res$centers gives K rows and D cols, where D is the dimension of the space
}


# # # # # # # # # 
# Regular Kmeans
# # # # # # # # #


regular.clustering <- function(X, K){
  kres <- cluster.kmeans(X, K)
  centers <- kres$centers
  fn <- function(x){
    res <- NULL
    for(i in seq(nrow(centers))){
      res <- c(res, sum((centers[i,]-x)^2))
    }
    return(as.vector(which.min(res)))
  }
  return(list(indicator=fn, centers=centers))
}


# # # # # # # # # # # # # # 
# Kmeans using projections
# # # # # # # # # # # # # #


# Projects the points on a random subspace of dimension d,
# and returns a function cfn that takes a point in the original
# scale, and cluster it. 
# Args:
# d is dimension of projected space
# K is number of clusters
# returns: an indicator function for which cluster
#          a certain point is in, and the centers (in
#          in the projected space)
.proj.clustering <- function(X, d, K){
  P <- get.subspace.D(ncol(X), d)
  pX <- t(P %*% t(X))
  kres = cluster.kmeans(pX, K)
  centers=kres$centers
  fn <- function(x){
    px <- as.vector(P %*% x)
    res <- NULL
    for(i in seq(1, nrow(centers))){
      res <- c(res, sum((centers[i,] - px)^2))
    }
    return(as.vector(which.min(res)))
  }
  return(list(cfn=fn, centers=centers))
}



# Computes the center of the clusters in the original scale
# Args:
# X is a matrix of data (each column is a data point)
# cfn an indicator function for the clusters
# K the number of clusters
get.augmented.centers <- function(X, cfn, K){
  clusters <- vector("list", K)

  for(i in seq(nrow(X))){
    clust.idx <- cfn(X[i,])
    clusters[[clust.idx]] <- rbind(clusters[[clust.idx]], X[i,])
  }

  centers <- NULL
  for(j in seq(K)){
    centers <- rbind(centers, colMeans(clusters[[j]]))
  }
  
  return(centers)
}


proj.clustering <- function(X, d, K){
  cfn <- .proj.clustering(X, d, K)
  centers <- get.augmented.centers(X, cfn, K)
  return(list(indicator=cfn, centers=centers))
}


# # # # # # # # # # # # # # 
# Density based-clustering
# # # # # # # # # # # # # #

detect.size.drop <- function(ds){
  n.clust <- max(unique(ds$cluster))
  size.clust <- NULL
  for(i in seq(n.clust)){
    size.clust <- c(size.clust, sum(ds$cluster==i))
  }
  sorted.clust <- sort(size.clust, decreasing=T, index.return=T)
  drops <- NULL
  for(i in seq(2,n.clust)){
    drops[i] <- sorted.clust$x[i] / sorted.clust$x[i-1]
  }
  print(drops)
  min.drop <- which.min(drops)
  high.clust <- sorted.clust$ix[seq(1,min.drop-1)]
  low.clust <- sorted.clust$ix[seq(min.drop, n.clust)]
  return(list(high=high.clust, low=low.clust))
}


thin.high.density <- function(draws, ds, high, thin.p=0.3, remove.noise=T){
  idx.0 <- which(ds$cluster==0)

  exclude.idx <- vector('list', length(high))

  for(j in seq(length(high))){
    exclude.idx[[j]] <- which(ds$cluster == high[j])
  }

  if(remove.noise){
    tmp.draws <- draws[-c(idx.0,do.call('c', exclude.idx)),]
  } else {
    tmp.draws <- draws[-do.call('c', exclude.idx),]
  }
  
  n <- round(thin.p*length(exclude.idx[[length(high)]]))

  print(n)

  for(j in seq(length(high))){
    tmp.draws <- rbind(tmp.draws, draws[sample(exclude.idx[[j]], n),])
  }
  return(tmp.draws)
}


smart.thinning <- function(draws, thin.p=0.2, scan.eps=0.3){
  niter <- ceiling(draws/7000)
  shuffled.idx <- sample(seq(nrow(draws)), nrow(draws))
  shuffled.draws <- draws[shuffled.idx,]
  res.draws <- NULL

  breakpoints <- seq(1, nrow(shuffled.draws), 7000)

  print(paste("total number of iterations:", length(breakpoints)-1))
  
  for(i in seq(length(breakpoints)-1)){
    print(paste(".. iteration:", i, "| length of slice:", breakpoints[i+1] - breakpoints[i]))
    tmp.draws <- shuffled.draws[seq(breakpoints[i], breakpoints[i+1]),]
    ds <- dbscan(tmp.draws, scan.eps)
    high <- detect.size.drop(ds)$high
    thin.draws <- thin.high.density(tmp.draws, ds, high, thin.p, remove.noise=F)
    res.draws <- rbind(res.draws, thin.draws)
  }
  return(res.draws)
}


thin.high.density.old <- function(draws, ds, thin.p=0.3){
  idx.0 <- which(ds$cluster==0)
  idx.1 <- which(ds$cluster==1)
  idx.2 <- which(ds$cluster==2)
  idx.3 <- which(ds$cluster==3)

  tmp.draws <- draws[-c(idx.0,idx.1,idx.2,idx.3),]

  #n.1 <- round(thin.p*length(idx.1))
  #n.2 <- round(thin.p*length(idx.2))
  n.3 <- round(thin.p*length(idx.3))

  #print(n.1)
  #print(n.2)
  print(n.3)
  
  tmp.draws <- rbind(tmp.draws, draws[sample(idx.1, n.3),])
  tmp.draws <- rbind(tmp.draws, draws[sample(idx.2, n.3),])
  tmp.draws <- rbind(tmp.draws, draws[sample(idx.3, n.3),])

  return(tmp.draws)
}


# # # #
# Spectral clustering stuff
# # # #

# this one is ~ 120 times faster than get.K.hat.slow
get.K.hat <- function(X, Lvec, nsub=500){
  Y.idx <- sample(seq(1, nrow(X)), nsub)
  Y <- X[Y.idx,]
  all.rows <- lapply(seq(nsub), function(i) Lvec(Y, Y[i,]))
  return(list(K.hat=do.call('rbind', all.rows), Y=Y))
}

get.K.hat.slow <- function(X, L, nsub=500){
  Y.idx <- sample(seq(1, nrow(X)), nsub)
  Y <- X[Y.idx,]
  get.row <- function(i) sapply(seq(nsub), function(j) L(Y[j,], Y[i,]))
  all.rows <- lapply(seq(nsub), get.row)
  return(list(K.hat=do.call('rbind', all.rows), Y=Y))
}


.get.closeness.centers <- function(Y, K.hat, clusters, K){
  centers <- NULL
  for(i in seq(K)){
    clust.idx <- which(clusters==i)
    sub.K.hat <- K.hat[clust.idx, clust.idx]
    #sub.K.hat <- t(apply(sub.K.hat, 1, function(x) ifelse(x == 0, 1/min(x), 1/x)))
    #center.idx <- which.min(colSums(sub.K.hat))
    center.idx <- which.max(colSums(sub.K.hat))
    centers <- rbind(centers, Y[clust.idx[center.idx],])
  }
  return(centers)
}

.get.closeness.centers.2 <- function(Y, K.hat, clusters, K, n.dist.0=5){
  centers <- NULL
  for(i in seq(K)){
    clust.idx <- which(clusters==i)
    n.dist <- min(length(clust.idx), n.dist.0)
    sub.K.hat <- K.hat[clust.idx, clust.idx]
    sub.K.hat <- t(apply(sub.K.hat, 1, function(x) ifelse(x == 0, min(x), x)))
    center.idx <- which.max(apply(sub.K.hat, 2, function(x) sum(sort(x)[1:n.dist])))
    centers <- rbind(centers, Y[clust.idx[center.idx],])
  }
  return(centers)
}

.get.closeness.centers.3.closure <- function(target.fn){
    .get.closeness.centers.3 <- function(Y, K.hat, clusters, K){
        centers <- NULL
        for(i in seq(K)){
            clust.idx <- which(clusters==i)
            max.dens.idx <- which.max(apply(Y[clust.idx,], 1, target.fn))
            centers <- rbind(centers, Y[clust.idx[max.dens.idx],])
        }
        return(centers)
    }
    return(.get.closeness.centers.3)
}

.do.spec.clustering <- function(K.hat, K){
  D.sq.inv.vec <- 1/sqrt(rowSums(K.hat))
  if(any(is.na(D.sq.inv.vec))){
    D.sq.inv.vec[is.na(D.sq.inv.vec)] <- 0
    warning("Some entries where set to 0 in D^(-1/2)")
  }
  D.sq.inv <- diag(D.sq.inv.vec)
  K.L <- D.sq.inv %*% K.hat %*% D.sq.inv
  
  if(any(is.na(K.L)) || any(is.infinite(K.L))){
    browser()
  }
  
  eigen.res <- eigen(K.L)
  eigenvecs <- eigen.res$vectors[,seq(1,K)]
  eigenvals <- eigen.res$values[seq(1,K)]
  
  
  if(is.complex(eigenvecs)){
    eigenvecs <- matrix(as.double(eigenvecs), nrow=nrow(eigenvecs))
    warning("complex eigenvectors were converted to real (imaginary part discarded)")
  }
  
  n.eigenvecs <- t(apply(eigenvecs, 1, function(x) x/sqrt(sum(x^2))))
  #save(n.eigenvecs, K.hat, K, K.L, file="debugfile.rdata")
  kmeans.res <- kmeans(n.eigenvecs, K, iter.max=400, nstart=15)
  return(list(cluster=kmeans.res$cluster, centers=kmeans.res$centers, n.eigenvecs=n.eigenvecs, eigenvals=eigenvals))
}

.get.cluster.fn <- function(centers, L, Y, eigenvecs){
  
  fn <- function(x){
    x.dist <- L(Y,x)    
    if(any(is.na(x.dist))){
      x.dist[is.na(x.dist)] <- 0
      warning("some NA entries where set to 0")
    }
    norm.x.dist <- x.dist / sum(x.dist)
    px <- as.vector(norm.x.dist %*% eigenvecs)
    npx <- px / sqrt(sum(px^2))
    dist.centers <- sapply(seq(nrow(centers)), function(i) sum( (npx-centers[i,])^2 ))
    
    return(which.min(dist.centers))
  }
  return(fn)
}


# note it is the responsibility of the user to check that L, get.K.hat, and cluster.method  work well together
# both should be vectorized or not.

spectral.clustering <- function(X, K, L, get.K.hat,
                                center.method=.get.closeness.centers.2,
                                cluster.method=.get.cluster.fn, nsub=100){

  # getting K.hat and Y
  K.res <- get.K.hat(X, L, nsub)
  K.hat <- K.res$K.hat
  Y <- K.res$Y
  #save(Y, file="debugfile-1.rdata")
  
  # the actual spectral clustering
  clust.res <- .do.spec.clustering(K.hat, K)
  
  # creating the cluster indicator function
  fn <- cluster.method(clust.res$centers, L, Y, clust.res$n.eigenvecs)

  # compute the centers used to initialize next mcmc iteration
  centers <- center.method(Y, K.hat, clust.res$cluster, K)
  
  return(list(indicator=fn, centers=centers, eigenvalues=clust.res$eigenvals, n.eigenvecs=clust.res$n.eigenvecs, Y=Y))
}

