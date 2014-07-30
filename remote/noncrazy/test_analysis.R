set.seed(123, "L'Ecuyer")

require(mvtnorm)
require(parallel)
require(methods)
source("../../mh.R", chdir=T)
source("../../bridge-sampling.R", chdir=T)
load("results/exploration.rdata")

indicator.fn <- sc$indicator


# the regular sampling part

run.simul <- function(j){
    res <- run.mv.mh(sc$centers[1,], 10000, rprop.1, dprop.1, constrained.targets[[1]])
    res.1.means <- colMeans(res$X)
    res <- run.mv.mh(sc$centers[2,], 10000, rprop.1, dprop.1, constrained.targets[[2]])
    res.2.means <- colMeans(res$X)
    res <- run.mv.mh(sc$centers[3,], 10000, rprop.1, dprop.1, constrained.targets[[3]])
    res.3.means <- colMeans(res$X)
    return(ws[1] * res.1.means + ws[2] * res.2.means + ws[3] * res.3.means)
}


res <- do.call('rbind', mclapply(seq(32), run.simul, mc.cores=32))

save(res, file="results/simul.rdata")
