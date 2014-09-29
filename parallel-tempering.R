
parallel.tempering <- function(targets.ls, single.rproposals.ls, single.dproposals.ls, init.states.mat, niter){

    nchains <- length(targets.ls)
    
    X.ls <- vector('list', nchains)
    L.ls <- vector('list', nchains)
    R.ls <- vector('list', nchains)
    
    X.tm1 <- init.states.mat
    acceptance.count <- rep(0, nchains)
    swap.acceptance.count <- 0
    
    for(i in seq(niter)){
        # regular state update
        for(j in seq(nchains)){
            L.t <- single.rproposals.ls[[j]](X.tm1[j,])
            r.num <- targets.ls[[j]](L.t) * single.dproposals.ls[[j]](X.tm1[j,], L.t)
            r.denom <- targets.ls[[j]](X.tm1[j,]) * single.dproposals.ls[[j]](L.t, X.tm1[j,])
            r <- r.num/r.denom
            X.ls[[j]] <- rbind(X.ls[[j]], X.tm1[j,])
            L.ls[[j]] <- rbind(L.ls[[j]], L.t)
            R.ls[[j]] <- c(R.ls[[j]], min(r,1))
            
            if(runif(1) < r){
                acceptance.count[j] <- acceptance.count[j] + 1
                X.tm1[j,] <- L.t
            }
        }

        # state swapping update
        indx <- sample(seq(nchains), replace=F)
        r.num <- targets.ls[[indx[1]]](X.tm1[indx[2],]) * targets.ls[[indx[2]]](X.tm1[indx[1],])
        r.denom <- targets.ls[[indx[1]]](X.tm1[indx[1],]) * targets.ls[[indx[2]]](X.tm1[indx[2],])
        if(runif(1) < r.num / r.denom){
            swap.acceptance.count <- swap.acceptance.count + 1
            tmp <- X.tm1[indx[1],]
            X.tm1[indx[1],] <-  X.tm1[indx[2],]
            X.tm1[indx[2],] <- tmp
        }

        for(j in seq(nchains)){
            X.ls[[j]][i,] <- X.tm1[j,]
        }
        
    }
    
    return(list(X.ls, L.ls, R.ls, swap.acceptance.count/niter, acceptance.count/niter))
}


parallel.tempering.notrace <- function(targets.ls, single.rproposals.ls, single.dproposals.ls, init.states.mat, niter){

    nchains <- length(targets.ls)
    
    X <- NULL
    
    X.tm1 <- init.states.mat
    acceptance.count <- rep(0, nchains)
    swap.acceptance.count <- 0
    
    for(i in seq(niter)){
        # regular state update
        for(j in seq(nchains)){
            L.t <- single.rproposals.ls[[j]](X.tm1[j,])
            r.num <- targets.ls[[j]](L.t) * single.dproposals.ls[[j]](X.tm1[j,], L.t)
            r.denom <- targets.ls[[j]](X.tm1[j,]) * single.dproposals.ls[[j]](L.t, X.tm1[j,])
            r <- r.num/r.denom
            
            if(runif(1) < r){
                acceptance.count[j] <- acceptance.count[j] + 1
                X.tm1[j,] <- L.t
            }
        }

        # state swapping update
        indx <- sample(seq(nchains), replace=F)
        r.num <- targets.ls[[indx[1]]](X.tm1[indx[2],]) * targets.ls[[indx[2]]](X.tm1[indx[1],])
        r.denom <- targets.ls[[indx[1]]](X.tm1[indx[1],]) * targets.ls[[indx[2]]](X.tm1[indx[2],])
        if(runif(1) < r.num / r.denom){
            swap.acceptance.count <- swap.acceptance.count + 1
            tmp <- X.tm1[indx[1],]
            X.tm1[indx[1],] <-  X.tm1[indx[2],]
            X.tm1[indx[2],] <- tmp
        }

        X <- rbind(X, X.tm1[1,])
    }
    
    return(list(X, swap.acceptance.count/niter, acceptance.count/niter))
}






###### Legacy stuff



parallel.tempering.old <- function(targets.ls, single.rproposals.ls, single.dproposals.ls, init.states.mat, niter){

    nchains <- length(targets.ls)
    
    X.ls <- vector('list', nchains)
    L.ls <- vector('list', nchains)
    R.ls <- vector('list', nchains)
    
    X.tm1 <- init.states.mat
    acceptance.count <- rep(0, nchains)
    swap.acceptance.count <- 0
    
    for(i in seq(niter)){
        # regular state update
        for(j in seq(nchains)){
            L.t <- single.rproposals.ls[[j]](X.tm1[j,])
            r.num <- targets.ls[[j]](L.t) * single.dproposals.ls[[j]](X.tm1[j,], L.t)
            r.denom <- targets.ls[[j]](X.tm1[j,]) * single.dproposals.ls[[j]](L.t, X.tm1[j,])
            r <- r.num/r.denom
            X.ls[[j]] <- rbind(X.ls[[j]], X.tm1[j,])
            L.ls[[j]] <- rbind(L.ls[[j]], L.t)
            R.ls[[j]] <- c(R.ls[[j]], min(r,1))
            
            if(runif(1) < r){
                acceptance.count[j] <- acceptance.count[j] + 1
                X.tm1[j,] <- L.t
            }
        }

        # state swapping update
        indx <- sample(seq(nchains), replace=F)
        r.num <- targets.ls[[indx[1]]](X.ls[[indx[2]]][i,]) * targets.ls[[indx[2]]](X.ls[[indx[1]]][i,])
        r.denom <- targets.ls[[indx[1]]](X.ls[[indx[1]]][i,]) * targets.ls[[indx[2]]](X.ls[[indx[2]]][i,])
        if(runif(1) < r.num / r.denom){
            swap.acceptance.count <- swap.acceptance.count + 1
            tmp <- X.ls[[indx[1]]][i,]
            X.ls[[indx[1]]][i,] <-  X.ls[[indx[2]]][i,]
            X.ls[[indx[2]]][i,] <- tmp
            X.tm1[indx[1],] <- X.ls[[indx[1]]][i,]
            X.tm1[indx[2],] <- X.ls[[indx[2]]][i,]
        }
    }
    
    return(list(X.ls, L.ls, R.ls, swap.acceptance.count/niter, acceptance.count/niter))
}

