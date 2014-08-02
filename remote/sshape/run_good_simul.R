require(methods)
require(mvtnorm)
require(parallel)
source("../../mh.R", chdir=T)
source("../../clustering.R", chdir=T)

load("results/exploration-sc.rdata")


run.simul <- function(i){
    res.sc.1 <- run.mv.mh(sc$centers[1,], 15000, rprop.1, dprop.1, sc.constrained.targets[[1]])
    res.sc.2 <- run.mv.mh(sc$centers[2,], 15000, rprop.1, dprop.1, sc.constrained.targets[[2]])
    return(ws.sc[1] * colMeans(res.sc.1$X) + ws.sc[2] * colMeans(res.sc.2$X))
}

res <- do.call('cbind', mclapply(seq(64), run.simul, mc.cores = 32))

save(res, file="results/simul-sc.rdata")
