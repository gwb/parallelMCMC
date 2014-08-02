require(methods)
require(mvtnorm)
require(parallel)
source("../../mh.R", chdir=T)
source("../../clustering.R", chdir=T)

load("results/exploration-vc.rdata")


run.simul <- function(i){
    res.vc.1 <- run.mv.mh(vc.centers[[1]], 15000, rprop.1, dprop.1, vc.constrained.targets[[1]])
    res.vc.2 <- run.mv.mh(vc.centers[[2]], 15000, rprop.1, dprop.1, vc.constrained.targets[[2]])
    return(ws.vc[1] * colMeans(res.vc.1$X) + ws.vc[2] * colMeans(res.vc.2$X))
}

res <- do.call('cbind', mclapply(seq(64), run.simul, mc.cores = 32))

save(res, file="results/simul-vc.rdata")
