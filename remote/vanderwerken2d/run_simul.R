set.seed(123, "L'Ecuyer")
require(methods)
require(mvtnorm)
require(parallel)

source("../../clustering.R", chdir=T)
source("../../mh.R", chdir=T)

load("results/exploration.rdata")


run.simul <- function(i){
    res <- run.mv.mh(sc$centers[1,], 5000, rprop.1, dprop.1, constrained.targets[[1]])
    res.1.means <- colMeans(res$X)
    res <- run.mv.mh(sc$centers[2,], 5000, rprop.1, dprop.1, constrained.targets[[2]])
    res.2.means <- colMeans(res$X)
    res <- run.mv.mh(sc$centers[3,], 5000, rprop.1, dprop.1, constrained.targets[[3]])
    res.3.means <- colMeans(res$X)
    res <- run.mv.mh(sc$centers[4,], 5000, rprop.1, dprop.1, constrained.targets[[4]])
    res.4.means <- colMeans(res$X)
    return(ws[1] * res.1.means + ws[2] * res.2.means + ws[3] * res.3.means + ws[4] * res.4.means)
}

res <- do.call('rbind', mclapply(seq(64), run.simul, mc.cores=32))

save(res, file="results/good.rdata")
