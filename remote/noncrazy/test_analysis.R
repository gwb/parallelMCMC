set.seed(123, "L'Ecuyer")

require(mvtnorm)
require(parallel)
require(methods)
source("../../mh.R", chdir=T)
source("../../bridge-sampling.R", chdir=T)
load("results/exploration.rdata")

indicator.fn <- sc$indicator



# the bridge part
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
    x2 <- rq2(5000)
    res.i <- bridge.sampling.very.fast(constrained.targets[[i]], dq2, xi, x2)
    return(list(res.i))
}

# the regular sampling part

run.simul <- function(j){
    res.1 <- run.mv.mh(sc$centers[1,], 10000, rprop.1, dprop.1, constrained.targets[[1]])
    res.2 <- run.mv.mh(sc$centers[2,], 10000, rprop.1, dprop.1, constrained.targets[[2]])
    res.3 <- run.mv.mh(sc$centers[3,], 10000, rprop.1, dprop.1, constrained.targets[[3]])
    res.ls <- list(res.1$X, res.2$X, res.3$X)
    wis <- unlist(lapply(seq(3), function(i) compute.weight.i(i, res.ls[[i]])))
    wis <- wis / sum(wis)
    comps <- lapply(seq(3), function(i) wis[i] * colMeans(res.ls[[i]]))
    res <- colSums(do.call('rbind', comps))
    return(res)
}


res <- do.call('rbind', mclapply(seq(32), run.simul, mc.cores=32))

save(res, file="results/simul.rdata")
