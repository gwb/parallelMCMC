set.seed(123, "L'Ecuyer")
require(methods)
require(parallel)
require(mvtnorm)
source("../../clustering.R", chdir=T)
source("../../parallel-tempering.R", chdir=T)
source("../../mh.R", chdir=T)


# The distribution below is not crazy, but it's pretty nasty
rfn <- function(){
    rho.1 <- 0.95
    rho.2 <- -0.95
    rho.3 <- 0
    sig.1 <- matrix(c(1, sqrt(8)*rho.1, sqrt(8)*rho.1, 8), nrow=2, byrow=T)
    sig.2 <- matrix(c(1, sqrt(8)*rho.2, sqrt(8)*rho.2, 8), nrow=2, byrow=T)
    sig.3 <- matrix(c(0.05, sqrt(0.05*0.1)*rho.3, sqrt(0.05*0.1)*rho.3, 0.1), nrow=2, byrow=T)
    rfn.1 <- function(){rmvnorm(1, c(-1.5, 2.5), sig.1)}
    rfn.2 <- function(){rmvnorm(1, c(1.5, 2.5), sig.2)}
    rfn.3 <- function(){rmvnorm(1, c(0,0), sig.3)}
    rfns <- list(rfn.1, rfn.2, rfn.3)
    #rfns <- list(rfn.1, rfn.2)
    idx <- sample(c(1,2,3),1, prob = c(0.4, 0.4, 0.2))
    #idx <- sample(c(1,2),1)
    return(rfns[[idx]]())
}

dfn <- function(x){
    rho.1 <- 0.95
    rho.2 <- -0.95
    rho.3 <- 0
    sig.1 <- matrix(c(1, sqrt(8)*rho.1, sqrt(8)*rho.1, 8), nrow=2, byrow=T)
    sig.2 <- matrix(c(1, sqrt(8)*rho.2, sqrt(8)*rho.2, 8), nrow=2, byrow=T)
    sig.3 <- matrix(c(0.05, sqrt(0.05*0.1)*rho.3, sqrt(0.05*0.1)*rho.3, 0.1), nrow=2, byrow=T)

    comp.1 <- 0.4 * dmvnorm(x, c(-1.5, 2.5), sig.1)
    comp.2 <- 0.4 * dmvnorm(x, c(1.5, 2.5), sig.2)
    comp.3 <- 0.2 *dmvnorm(x, c(0, 0), sig.3)
    
    return(comp.1 + comp.2 + comp.3)
}


# Exploration with parallel tempering
dfn.1 <- function(x) { dfn(x) }
dfn.2 <- function(x) { dfn(x)^(0.7) }
dfn.3 <- function(x) { dfn(x)^(0.4) }
dfn.4 <- function(x) { dfn(x)^(0.25) }


rprop.1 <- function(ostate){
    return(rmvnorm(1, ostate, 0.05 * diag(2)))
}
rprop.2 <- function(ostate){
    return(rmvnorm(1, ostate, 0.1 * diag(2)))
}
rprop.3 <- function(ostate){
    return(rmvnorm(1, ostate, 0.3 * diag(2)))
}
rprop.4 <- function(ostate){
    return(rmvnorm(1, ostate, 0.5 * diag(2)))
}

dprop.1 <- function(nstate, ostate){
    return(dmvnorm(nstate, ostate, 0.05 * diag(2)))
}
dprop.2 <- function(nstate, ostate){
    return(dmvnorm(nstate, ostate, 0.1 * diag(2)))
}
dprop.3 <- function(nstate, ostate){
    return(dmvnorm(nstate, ostate, 0.3 * diag(2)))
}
dprop.4 <- function(nstate, ostate){
    return(dmvnorm(nstate, ostate, 0.5 * diag(2)))
}



# preparing parameters for parallel tempering
targets.ls <- list(dfn.1, dfn.2, dfn.3, dfn.4)
single.rproposals.ls <- list(rprop.1, rprop.2, rprop.3, rprop.4)
single.dproposals.ls <- list(dprop.1, dprop.2, dprop.3, dprop.4)
init.states.mat <- matrix(rep(0,8), nrow=4)


run.simul <- function(j){
    res <- parallel.tempering(targets.ls, single.rproposals.ls, single.dproposals.ls, init.states.mat, 10000)
    states <- res[[1]][[1]]
    return(colMeans(states))
}


res <- do.call('rbind', mclapply(seq(32), run.simul, mc.cores=32))

save(res, file="results/parallel_tempering.rdata")
