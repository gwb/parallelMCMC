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
    res.1 <- run.mv.mh(sc$centers[1,], 10000, rprop.1, dprop.1, constrained.targets[[1]])
    res.2 <- run.mv.mh(sc$centers[2,], 10000, rprop.1, dprop.1, constrained.targets[[2]])
    res.3 <- run.mv.mh(sc$centers[3,], 10000, rprop.1, dprop.1, constrained.targets[[3]])
    return(ws[1] * colMeans(res.1$X) + ws[2] * colMeans(res.2$X) + ws[3] * colMeans(res.3$X))
}


res <- do.call('rbind', mclapply(seq(32), run.simul, mc.cores=32))

save(res, file="results/simul.rdata")
