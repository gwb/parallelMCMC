source('../../clustering.R', chdir=T)


E <- function(M){
    res <- 0
    for(i in seq(0,N-1)){
        for(j in seq(0,N-1)){
            idx.down <- 1+c( (i+1) %% N, j)
            idx.right <- 1+c(i, (j+1) %% N)
            res <- res + J * M[t(1+c(i,j))] * M[t(idx.right)]
            res <- res + J * M[t(1+c(i,j))] * M[t(idx.down)]
        }
    }
    return(-res)
}

unnorm.dens.vec <- function(M, Beta=1/(T*k.B), N.spin=N){
    return(exp(- E(matrix(M, nrow=N.spin)) * Beta))
}

get.dist <- function(M1, M2){
    #delta.M <- min(floor(N^2/2), N^2 - sum(M1 == M2))
    delta.M <- N^2 -sum(M1==M2)
    return(exp(-delta.M^2))
}

Q <- function(M.vec.1, M.vec.2, size=N){
    if(!is.null(nrow(M.vec.1))){
        tmp <- M.vec.2
        M.vec.2 <- M.vec.1
        M.vec.1 <- tmp
    }
    if(is.null(nrow(M.vec.2))){
        res <- get.dist(matrix(M.vec.2, nrow=size), matrix(M.vec.1, nrow=size))
    } else {
        res <- NULL
        for(i in seq(nrow(M.vec.2))){
            res <- c(res, get.dist(matrix(M.vec.1, nrow=size),
                                   matrix(M.vec.2[i,], nrow=size)))
        }
    }

    return(res)
}

D.E <- function(M, idx){
    tmp.idx <- idx - 1
    idx.down <- 1 + c( (tmp.idx[1]+1) %% N, tmp.idx[2] )
    idx.up <- 1 + c( (tmp.idx[1]-1) %% N, tmp.idx[2] )
    idx.right <- 1 + c( tmp.idx[1] , (tmp.idx[2]+1) %% N)
    idx.left <- 1 + c( tmp.idx[1], (tmp.idx[2]-1) %% N)
    return(2 * J * M[t(idx)] * sum(M[rbind(idx.down, idx.up, idx.right, idx.left)]))
}

prop.spin <- function(N){
    n <- sample(seq(0, N^2 - 1), 1)
    i <- floor(n/N) + 1
    j <- n %% N + 1
    return(c(i,j))
}

run.mh.modif <- function(M.0, prop.spin, get.A.i, Beta=1/(T*k.B), niter=100, N.spin=N){
    samples.ls <- vector('list', niter+1)
    samples.ls[[1]] <- M.0
    M.TM1 <- M.0
    n.accepted <- 0
    
    for(i in seq(niter)){
        idx.prop <- prop.spin(N)
        #D.E.prop <- D.E(M.TM1, idx.prop)
        #A.i <- min(1, exp(- Beta * D.E.prop))
        A.i <- get.A.i(M.TM1, idx.prop, Beta)
        if(runif(1) <= A.i){
            M.TM1[t(idx.prop)] <- -M.TM1[t(idx.prop)]
            n.accepted <- n.accepted + 1
        }
        samples.ls[[i+1]] <- M.TM1
    }

    return(list(samples.ls, n.accepted / niter))
}


run.mh.light.modif <- function(M.0, prop.spin, get.A.i, Beta=1/(T*k.B), niter=100, N.spin=N){
    magnet.ls <- vector('numeric', niter+1)
    magnet.ls[1] <- sum(M.0) / N.spin^2
    M.TM1 <- M.0
    n.accepted <- 0
    
    for(i in seq(niter)){
        idx.prop <- prop.spin(N)
        #D.E.prop <- D.E(M.TM1, idx.prop)
        #A.i <- min(1, exp(- Beta * D.E.prop))
        A.i <- get.A.i(M.TM1, idx.prop, Beta)
        if(runif(1) <= A.i){
            M.TM1[t(idx.prop)] <- -M.TM1[t(idx.prop)]
            n.accepted <- n.accepted + 1
        }
        magnet.ls[i+1] <- sum(M.TM1) / N.spin^2
    }

    return(list(magnet.ls, n.accepted / niter))
}

get.magnetization.distr <- function(samples.ls){
    magnet.dico <- cbind(seq(-1,1, by=2/N^2),0)
    for(i in seq(length(samples.ls))){
        magnet <- sum(samples.ls[[i]]) / N^2
        magnet.idx <- which(abs(magnet.dico[,1]-magnet) < 10^(-8))
        magnet.dico[magnet.idx,2] <- magnet.dico[magnet.idx,2] + 1
    }
    return(magnet.dico)
}

get.magnetization.distr.light <- function(magnet.ls){
    magnet.dico <- cbind(seq(-1,1, by=2/N^2),0)
    for(i in seq(length(magnet.ls))){
        magnet.idx <- which(abs(magnet.dico[,1]-magnet.ls[i]) < 10^(-8))
        magnet.dico[magnet.idx,2] <- magnet.dico[magnet.idx,2] + 1
    }
    return(magnet.dico)
}


.get.constrained.Ai <- function(constraint.fn, original.Ai){
    force(constraint.fn)
    constrained.Ai <- function(M, idx, bet){
        M.star <- M
        M.star[t(idx)] <- -M[t(idx)]
        cstr.fn.x <- constraint.fn(as.vector(M.star))
        if(length(cstr.fn.x) == 0){
            warning("something didn't go well")
            return(0)
        }
        if(cstr.fn.x){
            return(original.Ai(M, idx, bet))
        } else {
            return(0)
        }
    }
    return(constrained.Ai)
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
get.constrained.Ai <- function(K, cfn, original.Ai){
  res <- vector('list',K)
  for(i in seq(K)){
    res[[i]] <- .get.constrained.Ai(.get.constraint.fn(i, cfn), original.Ai)
    force(res)
  }
  return(res)
}



.get.A.i <- function(M, idx, bet){
    D.E.prop <- D.E(M, idx)
    A.i <- min(1, exp(- bet * D.E.prop))
    return(A.i)
}


### bridge sampling section

get.n.u.X <- function(res.mh, cutoff=0.99){
    X <- do.call('rbind', lapply(res.mh[[1]], function(x) as.vector(x)))
    u.X <- unique(X)
    res.pm <- NULL
    for(i in seq(nrow(u.X))){
        res.pm <- c(res.pm,
                    sum(apply(X, 1, function(x) all(x == u.X[i,]))))
    }

    res.pm <- res.pm / sum(res.pm)
    res.pm.sorted <- sort(res.pm, decreasing=TRUE, index.return=T)    
    idx <- which(cumsum(res.pm.sorted$x) >= cutoff)[1]
    keep.idx <- res.pm.sorted$ix[seq(idx)]
    keep.pm <- res.pm.sorted$x[seq(idx)]
    n.u.X <- u.X[keep.idx,]
    #browser()
    return(list(n.u.X, keep.pm))
}

get.r.prop <- function(n.u.X, keep.pm){
    r.prop <- function(){
        idx <- sample(seq(nrow(n.u.X)), 1, prob=keep.pm/sum(keep.pm))
        return(n.u.X[idx,])
    }
    return(r.prop)
}

get.d.prop <- function(n.u.X, keep.pm){
    d.prop <- function(x.prop){
        lookup <- apply(n.u.X, 1, function(x) all(x == x.prop))
        if(!any(lookup)){
            return(0)
        } else {
            idx <- which(lookup)
            return(keep.pm[idx])
        }
    }
    return(d.prop)
}


get.proposals <- function(res.mh, cutoff=0.99){
    res.nuX <- get.n.u.X(res.mh)
    n.u.X <- res.nuX[[1]]
    keep.pm <- res.nuX[[2]]

    r.prop <- get.r.prop(n.u.X, keep.pm)
    d.prop <- get.d.prop(n.u.X, keep.pm)

    return(list(r.prop, d.prop))
}


get.swap.A.i <- function(M.i, Beta.i, M.j, Beta.j){
    return( min(1, exp( (Beta.i - Beta.j) * (E(M.i) - E(M.j)))) )
}

parallel.tempering <- function(M.0.ls, prop.spin, get.A.i, get.swap.A.i, Beta.ls, niter=100, N.spin=N){


    nchains <- length(M.0.ls)
    samples.ls <- vector('list', nchains)
    #M.TM1 <- do.call('rbind', lapply(M.0.ls, function(x) as.vector(x)))
    M.TM1.ls <- M.0.ls
    n.accepted.ls <- vector('numeric', nchains)
    n.swap <- 0

    for(n in seq(niter)){
        ## Regular State Updates
        for(i in seq(nchains)){
            #if(i == 4){browser()}
            idx.prop <- prop.spin(N.spin)
            A.i <- get.A.i(M.TM1.ls[[i]], idx.prop, Beta.ls[i])
            if(runif(1) <= A.i){
                M.TM1.ls[[i]][t(idx.prop)] <- -M.TM1.ls[[i]][t(idx.prop)]
                n.accepted.ls[i] <- n.accepted.ls[i] + 1
            }
            samples.ls[[i]] <- rbind(samples.ls[[i]], as.vector(M.TM1.ls[[i]])) 
        }
        
        ## State Swapping Updates
        indx <- sample(seq(nchains), replace=F)
        if(runif(1) <= get.swap.A.i(M.TM1.ls[[indx[1]]], Beta.ls[indx[1]], M.TM1.ls[[indx[2]]], Beta.ls[indx[2]])){
            
            n.swap <- n.swap + 1
            tmp <- M.TM1.ls[[indx[1]]]
            M.TM1.ls[[indx[1]]] <- M.TM1.ls[[indx[2]]]
            M.TM1.ls[[indx[2]]] <- tmp
        }
        
        for(j in seq(nchains)){
            samples.ls[[j]][n,] <- as.vector(M.TM1.ls[[j]])
        }
    }
    return(list(samples.ls, n.accepted.ls / niter, n.swap / niter))
}
