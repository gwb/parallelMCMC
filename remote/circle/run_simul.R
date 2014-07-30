set.seed(123, "L'Ecuyer")
require(methods)
require(mvtnorm)
require(parallel)

source("../../clustering.R", chdir=T)
source("../../mh.R", chdir=T)

load("results/exploration.rdata")


run.simul <- function(i){
    res.1 <- run.mv.mh(sc$centers[1,], 10000, rprop.1, dprop.1, constrained.targets[[1]])
    res.2 <- run.mv.mh(sc$centers[2,], 10000, rprop.1, dprop.1, constrained.targets[[2]])
    return(ws[1] * colMeans(res.1$X) + ws[2] * colMeans(res.2$X))
}

res <- do.call('rbind', mclapply(seq(64), run.simul, mc.cores=32))

save(res, file="results/good.rdata")
