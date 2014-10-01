#
# The algorithm for spectral partitioning has been rewritten to reflect two major changes:
# 1) the use of the proposal kernel Q rather than the transition kernel L as a distance metric
# 2) the use of Bengio et al's method for out of sample extension of clustering
#

cluster.kmeans <- function(X,K,iter.max=50,nstart=20){
  # each row of X should be a point. so the number of columns of X should be the dimension of the space
  res <- kmeans(X,K,iter.max=iter.max,nstart=nstart)
  return(res) #res$centers gives K rows and D cols, where D is the dimension of the space
}




# # # # # # # # # # # # # # # # # # # # #
# Spectral Partitioning | Slow version  # 
# # # # # # # # # # # # # # # # # # # # #

# obtaining the transformed kernel
# (equivalent to first computing the Gram matrix, and then normalizing with t(D^-1/2) M D^-1/2

get.transformed.kernel <- function(X, Q){
    K <- function(new.state, old.state){
        av.1 <- sum(Q(new.state, X))/nrow(X)
        av.2 <- sum(Q(X, old.state))/nrow(X)
        return(Q(new.state, old.state)/sqrt(av.1*av.2))
    }
    return(K)
}


get.gram.matrix <- function(X, transformed.kernel){
    N <- nrow(X)
    M <- matrix(0, nrow=N, ncol=N)
    for(i in seq(2, N)){
        for(j in seq(1, i-1)){
            M[i,j] <- transformed.kernel(X[i,], X[j,])
            M[j,i] <- M[i,j]
        }
    }
    for(i in seq(N)){
        M[i,i] <- transformed.kernel(X[i,], X[i,])
    }
    return(M)
}

get.fk <- function(X, V, lambdas, transformed.kernel){
    n <- nrow(X)
    fk <- function(x, k){
        return( (sqrt(n) / lambdas[k]) * sum(sapply(seq(n), function(i) V[i,k] * transformed.kernel(x, X[i,]))))
    }
    return(fk)
}

get.embedding.fn <- function(X, V, lambdas, transformed.kernel){
    fk <- get.fk(X, V, lambdas, transformed.kernel)
    embedding.fn <- function(x){
        num.dim <- dim(V)[2]
        n <- nrow(X)
        return( (1/sqrt(n)) * sapply(seq(num.dim), function(i) fk(x, i)))
    }
    return(embedding.fn)
}


get.cluster.fn <- function(centers, embedding.fn){
    cluster.fn <- function(x){
        embedded.x <- embedding.fn(x)
        dist.centers <- sapply(seq(nrow(centers)), function(i) sum( (embedded.x - centers[i,])^2 ))
        return(which.min(dist.centers))
    }
    return(cluster.fn)
}

get.centers <- function(Y, clusters, nclust){
    res <- NULL
    for(i in seq(nclust)){
        idx.ls <- which(clusters == i)
        idx <- sample(idx.ls, 1)
        res <- rbind(res, Y[idx,])
    }
    return(res)
}

.do.spectral.clustering <- function(K.hat, nclust){

    eigen.res <- eigen(K.hat)
    eigenvecs <- eigen.res$vectors[,seq(1,nclust)]
    eigenvals <- eigen.res$values[seq(1,nclust)]

    if(is.complex(eigenvecs)){
        eigenvecs <- matrix(as.double(eigenvecs), nrow=nrow(eigenvecs))
        warning("complex eigenvectors were converted to real (imaginary part discarded)")
    }

    n.eigenvecs <- t(apply(eigenvecs, 1, function(x) x/sqrt(sum(x^2))))
    kmeans.res <- kmeans(n.eigenvecs, nclust, iter.max=1000, nstart=15)
    return(list(cluster=kmeans.res$cluster, centers=kmeans.res$centers, n.eigenvecs=n.eigenvecs, eigenvecs=eigenvecs, eigenvals=eigenvals))
}


spectral.clustering <- function(X, nclust, Q, nsub=100){

    # subsample X
    Y.idx <- sample(seq(1, nrow(X)), nsub)
    Y <- X[Y.idx,]

    # transform the kernel and obtain the gram matrix
    transformed.kernel <- get.transformed.kernel(Y, Q)
    K.hat <- get.gram.matrix(Y, transformed.kernel)
    
    # the actual spectral clustering
    clust.res <- .do.spectral.clustering(K.hat, nclust)

    # creating the cluster indicator function
    embedding.fn <- get.embedding.fn(Y, clust.res$eigenvecs, clust.res$eigenvals, transformed.kernel)
    cluster.indicator.fn <- get.cluster.fn(clust.res$centers, embedding.fn)

    # compute the centers used to initialize next mcmc iteration
    centers <- get.centers(Y, clust.res$cluster, nclust)

    return(list(indicator=cluster.indicator.fn, centers = centers, Y=Y, cluster = clust.res$cluster, embedding.fn=embedding.fn))
}



# # # # # # # # # # # # # # # # # # # # #
# Spectral Partitioning | Fast version  # 
# # # # # # # # # # # # # # # # # # # # #

# Note that the following has been optimized to work for symmetric kernels

get.vec.transformed.kernel <- function(X, Q){
    n <- nrow(X)
    av2 <- NULL
    for(i in seq(n)){
        av2 <- c(av2, sum(Q(X,X[i,])) / n)
    }
    
    K <- function(new.state, old.state=NULL){
        av.1 <- sum(Q(new.state, X)) / n
        if(!is.null(old.state)){
            av2.i <- sum(Q(X,old.state)) / n
            return(Q(old.state, new.state)/sqrt(av.1*av2.i))
        }
        return(Q(X, new.state) / sqrt(av.1*av2))
    }
    return(K)
}

get.vec.transformed.kernel.old <- function(X, Q){
    K <- function(new.state, old.state){
        av.1 <- sum(Q(new.state,X)) / nrow(X)
        if(!is.null(dim(old.state))){
            av.2 <- sapply(seq(nrow(old.state)), function(i) sum(Q(X, old.state[i,]))/nrow(X))
        } else {
            av.2 <- sum(Q(X, old.state))/nrow(X)
        }
        return(Q(old.state, new.state)/sqrt(av.1*av.2))
    }
    return(K)
}

get.fast.embedding.fn <- function(X, V, lambdas, transformed.kernel.vec){
    num.dim <- dim(V)[2]
    n <- nrow(X)
    embedding.fn <- function(x){
        KxXi <- transformed.kernel.vec(x)
        return( sapply(seq(num.dim), function(k) (1/lambdas[k]) * t(V[,k]) %*% KxXi) )
    }
    return(embedding.fn)
}


get.fast.embedding.fn.old <- function(X, V, lambdas, transformed.kernel.vec){
    num.dim <- dim(V)[2]
    n <- nrow(X)
    embedding.fn <- function(x){
        KxXi <- transformed.kernel.vec(x, X)
        return( sapply(seq(num.dim), function(k) (1/lambdas[k]) * t(V[,k]) %*% KxXi) )
    }
    return(embedding.fn)
}

spectral.clustering.vec <- function(X, nclust, Q, nsub=100){

    # subsample X
    Y.idx <- sample(seq(1, nrow(X)), nsub)
    Y <- X[Y.idx,]

    # transform the kernel and obtain the gram matrix
    transformed.kernel.vec <- get.vec.transformed.kernel(Y, Q)
    K.hat <- get.gram.matrix(Y, transformed.kernel.vec)
    
    # the actual spectral clustering
    clust.res <- .do.spectral.clustering(K.hat, nclust)

    # creating the cluster indicator function
    embedding.fn <- get.fast.embedding.fn(Y, clust.res$eigenvecs, clust.res$eigenvals, transformed.kernel.vec)
    cluster.indicator.fn <- get.cluster.fn(clust.res$centers, embedding.fn)

    # compute the centers used to initialize next mcmc iteration
    centers <- get.centers(Y, clust.res$cluster, nclust)

    return(list(indicator=cluster.indicator.fn, centers = centers, Y=Y, cluster = clust.res$cluster, embedding.fn=embedding.fn))
}

spectral.clustering.vec.old <- function(X, nclust, Q, nsub=100){

    # subsample X
    Y.idx <- sample(seq(1, nrow(X)), nsub)
    Y <- X[Y.idx,]

    # transform the kernel and obtain the gram matrix
    transformed.kernel.vec <- get.vec.transformed.kernel.old(Y, Q)
    K.hat <- get.gram.matrix(Y, transformed.kernel.vec)
    
    # the actual spectral clustering
    clust.res <- .do.spectral.clustering(K.hat, nclust)

    # creating the cluster indicator function
    embedding.fn <- get.fast.embedding.fn.old(Y, clust.res$eigenvecs, clust.res$eigenvals, transformed.kernel.vec)
    cluster.indicator.fn <- get.cluster.fn(clust.res$centers, embedding.fn)

    # compute the centers used to initialize next mcmc iteration
    centers <- get.centers(Y, clust.res$cluster, nclust)

    return(list(indicator=cluster.indicator.fn, centers = centers, Y=Y, cluster = clust.res$cluster, embedding.fn=embedding.fn))
}


# # # # # #
# Voronoi #
# # # # # # 

voronoi.clustering <- function(X, K){
  kres <- cluster.kmeans(X, K)
  centers <- kres$centers
  fn <- function(x){
    res <- NULL
    for(i in seq(nrow(centers))){
      res <- c(res, sum((centers[i,]-x)^2))
    }
    return(as.vector(which.min(res)))
  }
  return(list(indicator=fn, centers=centers, cluster=kres$cluster))
}
