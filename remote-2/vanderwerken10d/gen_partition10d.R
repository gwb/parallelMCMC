require(methods)
require(mvtnorm)
source("../../parallel-tempering.R", chdir=T)
source("../../mh.R", chdir=T)
source("../../clustering.R", chdir=T)
source("../../bridge-sampling.R", chdir=T)
source("../../graphics.R", chdir=T)
require(mvtnorm)

set.seed(123)


# Define target and all

C <- 4
d <- 10

sig.ls <- vector('list', C)
u.ls <- vector('list', C)

#w.ls <- c(0.1, 0.2, 0.3, 0.4)

sig.prop <- 0.2*diag(d)

for(i in seq(C)){
    sig.ls[[i]] <- 0.5*diag(d)
    u.ls[[i]] <- runif(d, min=-10, max=10)
}


# targets

dtarget <- function(x){
    w.ls <- c(0.1, 0.2, 0.3, 0.4)
      res <- 0
        for(i in seq(C)){
                res <- res + w.ls[i] * dmvnorm(x, u.ls[[i]], sig.ls[[i]])
            }
        return(res)
  }

rtarget <- function(n){
    w.ls <- c(0.1, 0.2, 0.3, 0.4)
      res <- NULL
        for(j in seq(n)){
                i <- sample(seq(C), 1, prob=w.ls)
                res <- rbind(res, rmvnorm(1, u.ls[[i]], sig.ls[[i]]))
            }
        return(res)
  }


dfn.1 <- function(x) {dtarget(x)}
dfn.2 <- function(x) {dtarget(x)^(0.75)}
dfn.3 <- function(x) {dtarget(x)^(0.5)}
dfn.4 <- function(x) {dtarget(x)^(0.1)}
dfn.5 <- function(x) {dtarget(x)^(0.05)}
dfn.6 <- function(x) {dtarget(x)^(0.01)}


dprop.1 <- function(nstate, ostate){
    return(dmvnorm(nstate, mean=ostate, sigma=sig.prop))
}
dprop.2 <- function(nstate, ostate){
        return(dmvnorm(nstate, ostate, 2 * sig.prop))
    }
dprop.3 <- function(nstate, ostate){
        return(dmvnorm(nstate, ostate, 5 * sig.prop))
    }
dprop.4 <- function(nstate, ostate){
        return(dmvnorm(nstate, ostate, 7 * sig.prop))
    }
dprop.5 <- function(nstate, ostate){
        return(dmvnorm(nstate, ostate, 10 * sig.prop))
    }
dprop.6 <- function(nstate, ostate){
        return(dmvnorm(nstate, ostate, 14 * sig.prop))
    }



rprop.1 <- function(ostate){
    return(rmvnorm(1, mean=ostate, sigma=sig.prop))
}
rprop.2 <- function(ostate){
        return(rmvnorm(1, ostate, 2 * sig.prop))
    }
rprop.3 <- function(ostate){
        return(rmvnorm(1, ostate, 5 * sig.prop))
    }
rprop.4 <- function(ostate){
        return(rmvnorm(1, ostate, 7 * sig.prop))
    }
rprop.5 <- function(ostate){
        return(rmvnorm(1, ostate, 10 * sig.prop))
    }
rprop.6 <- function(ostate){
        return(rmvnorm(1, ostate, 14 * sig.prop))
    }


dprop.1.vec <- function(nstate, ostate){
    if(!is.null(dim(ostate))){
        if(dim(ostate)[1] == 1){
            ostate <- as.vector(ostate)
        }
        return(dprop.1(ostate, nstate))
    } else if(!is.null(dim(nstate))){
        if(dim(nstate)[1] == 1){
            nstate <- as.vector(nstate)
        }
        return(dprop.1(nstate, ostate))
    } else {
        return(dprop.1(nstate, ostate))
    }
}



# Exploring with parallel tempering
targets.ls <- list(dfn.1, dfn.2, dfn.3, dfn.4, dfn.5, dfn.6)
single.rproposals.ls <- list(rprop.1, rprop.2, rprop.3, rprop.4, rprop.5, rprop.6)
single.dproposals.ls <- list(dprop.1, dprop.2, dprop.3, dprop.4, dprop.5, dprop.6)
init.states.mat <- matrix(rep(u.ls[[1]], 6), nrow=6, byrow=T)

set.seed(123)

#res <- parallel.tempering(targets.ls, single.rproposals.ls, single.dproposals.ls, init.states.mat, 20000)

res <- parallel.tempering.notrace(targets.ls, single.rproposals.ls, single.dproposals.ls, init.states.mat, 30000)
states <- res[[1]]

algo.dt <- plot.vhd.no.clusters(states)
ggplot(data=algo.dt[[1]], aes(x=x,y=y)) + geom_point() + facet_grid(.~proj)

# note: the error - 
# Error in do_one(nmeth) (from clustering.R#98) : NA/NaN/Inf in foreign function call (arg 1) can
# can be thrown if the variance of the proposal is too small

sc <- spectral.clustering.vec(states, 4, dprop.1.vec, nsub=500)

indicator.fn <- function(x) sc$indicator(as.vector(x))
constrained.targets <- get.constrained.targets(4, indicator.fn, dtarget)


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


save(C, d, ws, dtarget, dfn.1, rprop.1, dprop.1, dprop.1.vec, indicator.fn, constrained.targets, sc, u.ls, sig.ls, sig.prop, file="results/exploration.rdata")
