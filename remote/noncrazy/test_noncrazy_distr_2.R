require(mvtnorm)
source("../../clustering.R", chdir=T)
source("../../parallel-tempering.R", chdir=T)
source("../../mh.R", chdir=T)

set.seed(123)

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

# Visualizing the distribution
x <- do.call('rbind', lapply(seq(2000), function(i) rfn()))

#pdf(file="noncrazy.pdf")
#plot(x[,2] ~ x[,1])
#dev.off()


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


# Exploring with parallel tempering
targets.ls <- list(dfn.1, dfn.2, dfn.3, dfn.4)
single.rproposals.ls <- list(rprop.1, rprop.2, rprop.3, rprop.4)
single.dproposals.ls <- list(dprop.1, dprop.2, dprop.3, dprop.4)
init.states.mat <- matrix(rep(0,8), nrow=4)
res <- parallel.tempering(targets.ls, single.rproposals.ls, single.dproposals.ls, init.states.mat, 8000)


# Plotting what the chain of interest has explored
states <- res[[1]][[1]]
#plot(states[,2] ~ states[,1])


# Perform spectral clustering and get the cluster indicator function
get.K.hat.pruned.closure <- function(pprune){
    fn <- function(X, Lvec, nsub=500){
        return(get.K.hat.pruned(X, Lvec, nsub, pprune))
    }
    return(fn)
}

fn.get.K.hat.pruned <- get.K.hat.pruned.closure(0.04)

set.seed(456)

sc <- spectral.clustering(states, 3, dprop.1, fn.get.K.hat.pruned, nsub=700)
indicator.fn <- sc$indicator
#plot(sc$Y[,2] ~ sc$Y[,1], col=sc$clust.clust)

# Plots explored states with clusters in color
#clust.val <- apply(states, 1, indicator.fn)
#plot(states[,2] ~ states[,1], col=clust.val)


# running mh chains
constrained.targets <- get.constrained.targets(3, indicator.fn, dfn.1)

save(dfn, dfn.1, rprop.1, dprop.1, constrained.targets, sc, file="results/exploration.rdata")
