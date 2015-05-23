

# this function is the backbone for both the "exact" and "empirical"
# transition kernels. Q can be either the proposal or the transition
# kernel (leading to a proposal matrix, or transition matrix)
get.K <- function(X, Q, trim=TRUE){
    n <- nrow(X)
    K <- matrix(0, nrow=n, ncol=n)
    for(i in seq(1, n)){
        for(j in seq(1,n)){
            K[i,j] <- Q(X[i], X[j])
        }
    }
    K[is.na(K)] <- 0

    # normalize
    D.sqrt.inv <- diag(1/sqrt(rowSums(K)))
    K.norm <- D.sqrt.inv %*% K %*% D.sqrt.inv

    # Removes lines where all entries are 0
    # (these correspond to isolated points, and they mess things up,
    # by making spectral gap equal to zero, since they are unreachable)

    if(trim){
        idx.rm <- which(diag(K.norm) == 1)
        if(length(idx.rm)>0){
            K.norm <- K.norm[-idx.rm, -idx.rm]
        }
    }
        
    return(K.norm)
}

get.empirical.spectral.gaps <- function(X, mid, Q){
    Xv <- as.vector(X)
    X1 <- matrix(Xv[Xv < mid], ncol=1)
    X2 <- matrix(Xv[Xv > mid], ncol=1)
    K1 <- get.K(X1, Q)
    K2 <- get.K(X2, Q)
    res1 <- 1-eigen(K1, only.values=T)$values[2]
    res2 <- 1-eigen(K2, only.values=T)$values[2]
    if(any(is.complex(res1))){
        warning("complex eigenvalues where converted for K1")
        res1 <- as.double(res1)
    }
    if(any(is.complex(res2))){
        warning("complex eigenvalues where converted for K2")
        res2 <- as.double(res2)
    }
    return(c(res1, res2))
}

get.exact.spectral.gaps <- function(mu, mid, Q, n.discretize=300){
    X1 <- matrix(seq(-mu - 4*1, mid, length.out=n.discretize), ncol=1)
    X2 <- matrix(seq(mid, mu + 4*1, length.out=n.discretize), ncol=1)
    K1 <- get.K(X1, Q)
    K2 <- get.K(X2, Q)
    res1 <- 1-eigen(K1, only.values=T)$values[2]
    res2 <- 1-eigen(K2, only.values=T)$values[2]
    if(any(is.complex(res1))){
        warning("complex eigenvalues where converted for K1")
        res1 <- as.double(res1)
    }
    if(any(is.complex(res2))){
        warning("complex eigenvalues where converted for K2")
        res2 <- as.double(res2)
    }
    return(c(res1, res2))
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

get.partition.edge <- function(clusters, ys){
    y1 <- ys[clusters==1]
    y2 <- ys[clusters==2]
    t1 <- min(ys[clusters==1]) - max(ys[clusters==2])
    t2 <- min(ys[clusters==2]) - max(ys[clusters==1])
    
    if( t1 < 0 && t2 < 0){
        warning("overlapping clusters, special treatment required")
        ymax <- ifelse(max(y1) > max(y2), y1, y2)
        ymin <- ifelse(max(y1) > max(y2), y2, y1)
        ymax.overlap.mean <- mean(ymax[ymax < max(ymin)])
        ymin.overlap.mean <- mean(ymin[ymin > min(ymax)])
        if(is.na(ymax.overlap.mean) || is.na(ymin.overlap.mean)){
            warning("really bad partitioning, can't do anything")
            return(NA)
        }
        return(mean(ymax.overlap.mean, ymin.overlap.mean))
    } else{
        return(ifelse(t1>=0, t1, t2))
    }
}



get.midpoint <- function(n.mh, n.sub, n.burnin, rproposal, dproposal, dtarget){
    Xs <- run.mh(0, n.mh, rproposal, dproposal, dtarget, n.burnin)[[1]][,1]
    Ys <- matrix(sample(Xs, n.sub), ncol=1)
    K.hat <- get.K(Ys, dproposal, trim=F)
    clust.res <- .do.spectral.clustering(K.hat, 2)
    Ys.vec <- as.vector(Ys)

    mid <- get.partition.edge(clust.res[[1]], Ys.vec)
    return(list(Ys, mid, Xs, clust.res))
}
