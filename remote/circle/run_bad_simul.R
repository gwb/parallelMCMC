set.seed(123, "L'Ecuyer")
require(methods)
require(mvtnorm)
require(parallel)

source("../../mh.R", chdir=T)

load("results/exploration.rdata")


run.simul <- function(i){
    res.sim.1 <- run.mv.mh(sc$centers[1,], 10000, rprop.1, dprop.1, dcircle)
    res.sim.2 <- run.mv.mh(sc$centers[2,], 10000, rprop.1, dprop.1, dcircle)
    res.sim <- rbind(res.sim.1$X, res.sim.2$X)
    return(colMeans(res.sim))
}

res <- do.call('rbind', mclapply(seq(2), run.simul, mc.cores=2))

save(res, file="results/bad.rdata")
