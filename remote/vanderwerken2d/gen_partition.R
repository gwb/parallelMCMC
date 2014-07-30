require(methods)
require(mvtnorm)
source("../../parallel-tempering.R", chdir=T)
source("../../mh.R", chdir=T)
source("../../clustering.R", chdir=T)
source("../../bridge-sampling.R", chdir=T)

set.seed(123)

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

#x <- do.call('rbind', lapply(seq(1000), function(i) rvw()))
#plot(x[,2] ~ x[,1])


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
res <- parallel.tempering(targets.ls, single.rproposals.ls, single.dproposals.ls, init.states.mat, 8000)


states <- res[[1]][[1]]

#plot(states[,2]~states[,1])


sc <- spectral.clustering(states, 4, dprop.1, get.K.hat, nsub=700)
indicator.fn <- sc$indicator
plot(sc$Y[,2] ~ sc$Y[,1], col=sc$clust.clust)


# running mh chains
constrained.targets <- get.constrained.targets(4, indicator.fn, dfn.1)



get.r.centered.proposal <- function(u, sigmat){
    fn <- function(n){rmvt(n, delta=u, sigma=sigmat, df=10, type="shifted")}
    return(fn)
}

get.d.centered.proposal <- function(u,sigmat){
    fn <- function(x){dmvt(x, delta=u, sigma=sigmat, df=10, log=F)}
    return(fn)
}

compute.weight.i <- function(i, xi){
    u <- optim(sc$centers[i,], constrained.targets[[i]], control=list(fnscale=-1))$par
    sigmat <- cov(xi)
    rq2 <- get.r.centered.proposal(u, sigmat)
    dq2 <- get.d.centered.proposal(u, sigmat)
    x2 <- rq2(10000)
    res.i <- bridge.sampling.very.fast(constrained.targets[[i]], dq2, xi, x2)
    return(list(res.i))
}

run.simul <- function(){
    res.1 <- run.mv.mh(sc$centers[1,], 5000, rprop.1, dprop.1, constrained.targets[[1]])
    res.2 <- run.mv.mh(sc$centers[2,], 5000, rprop.1, dprop.1, constrained.targets[[2]])
    res.3 <- run.mv.mh(sc$centers[3,], 5000, rprop.1, dprop.1, constrained.targets[[3]])
    res.4 <- run.mv.mh(sc$centers[4,], 5000, rprop.1, dprop.1, constrained.targets[[4]])
    res.ls <- list(res.1$X, res.2$X, res.3$X, res.4$X)
    wis <- unlist(lapply(seq(4), function(i) compute.weight.i(i, res.ls[[i]])))
    wis <- wis / sum(wis)
    return(wis)
}


ws <- run.simul()

save(ws, dvw, dfn.1, rprop.1, dprop.1, constrained.targets, sc, u.ls, sig.ls, file="results/exploration.rdata")
