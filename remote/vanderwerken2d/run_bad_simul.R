set.seed(123, "L'Ecuyer")
require(methods)
require(mvtnorm)
require(parallel)
source("../../parallel-tempering.R", chdir=T)
source("../../mh.R", chdir=T)
source("../../clustering.R", chdir=T)
source("../../bridge-sampling.R", chdir=T)


u.ls <- list(c(3,3),
             c(7, -3),
             c(2,7),
             c(-5,0))

sig.ls <- list(matrix(c(1,0.2,0.2,1), nrow=2, byrow=T),
               matrix(c(2,-0.5,-0.5,0.5), nrow=2, byrow=T),
               matrix(c(1.3,0.3,0.3, 0.4), nrow=2, byrow=T),
               matrix(c(1,1,1,2.5), nrow=2, byrow=T))


rvw <- function(){
    idx <- sample(c(1,2,3,4), 1, prob=c(0.02, 0.2, 0.2, 0.58))
    return(rmvnorm(1, mean=u.ls[[idx]], sigma=sig.ls[[idx]]))
}

dvw <- function(x){
    comp.1 <- 0.02 * dmvnorm(x, mean=u.ls[[1]], sigma=sig.ls[[1]])
    comp.2 <- 0.2 * dmvnorm(x, mean=u.ls[[2]], sigma=sig.ls[[2]])
    comp.3 <- 0.2 * dmvnorm(x, mean=u.ls[[3]], sigma=sig.ls[[3]])
    comp.4 <- 0.58 * dmvnorm(x, mean=u.ls[[4]], sigma=sig.ls[[4]])
    return(comp.1 + comp.2 + comp.3 + comp.4)
}

dfn.1 <- function(x) { dvw(x) }
dfn.2 <- function(x) { dvw(x)^(0.7) }
dfn.3 <- function(x) { dvw(x)^(0.4) }
dfn.4 <- function(x) { dvw(x)^(0.25) }


rprop.1 <- function(ostate){
        return(rmvnorm(1, ostate, 0.8 * diag(2)))
    }
rprop.2 <- function(ostate){
        return(rmvnorm(1, ostate, 1 * diag(2)))
    }
rprop.3 <- function(ostate){
        return(rmvnorm(1, ostate, 2 * diag(2)))
    }
rprop.4 <- function(ostate){
        return(rmvnorm(1, ostate, 3 * diag(2)))
    }

dprop.1 <- function(nstate, ostate){
        return(dmvnorm(nstate, ostate, 0.8 * diag(2)))
    }
dprop.2 <- function(nstate, ostate){
        return(dmvnorm(nstate, ostate, 1 * diag(2)))
    }
dprop.3 <- function(nstate, ostate){
        return(dmvnorm(nstate, ostate, 2 * diag(2)))
    }
dprop.4 <- function(nstate, ostate){
        return(dmvnorm(nstate, ostate, 3 * diag(2)))
    }



# Exploring with parallel tempering
targets.ls <- list(dfn.1, dfn.2, dfn.3, dfn.4)
single.rproposals.ls <- list(rprop.1, rprop.2, rprop.3, rprop.4)
single.dproposals.ls <- list(dprop.1, dprop.2, dprop.3, dprop.4)
init.states.mat <- matrix(rep(c(-5,0),4), nrow=4, byrow=T)


run.simul <- function(i){
    states <- NULL
    for(i in seq(1,4)){
        res <- parallel.tempering(targets.ls, single.rproposals.ls, single.dproposals.ls, init.states.mat, 5000)
        states <- rbind(states, res[[1]][[1]])
    }
    return(colMeans(states))
}

res <- do.call('rbind', mclapply(seq(64), run.simul, mc.cores=32))

save(res, file="results/bad.rdata")
