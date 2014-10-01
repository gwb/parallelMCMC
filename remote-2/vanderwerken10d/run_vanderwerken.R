require(methods)
require(mvtnorm)
source("../../parallel-tempering.R", chdir=T)
source("../../mh.R", chdir=T)
source("../../clustering.R", chdir=T)
source("../../bridge-sampling.R", chdir=T)

args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[1])

load('results/exploration.rdata')



true.mean <- 0.1 * u.ls[[1]] + 0.2 * u.ls[[2]] + 0.3 * u.ls[[3]] + 0.4 * u.ls[[4]]

run.simul <- function(){
    res.1 <- run.mv.mh(sc$centers[1,], 7000, rprop.1, dprop.1, constrained.targets[[1]])
    res.2 <- run.mv.mh(sc$centers[2,], 7000, rprop.1, dprop.1, constrained.targets[[2]])
    res.3 <- run.mv.mh(sc$centers[3,], 7000, rprop.1, dprop.1, constrained.targets[[3]])
    res.4 <- run.mv.mh(sc$centers[4,], 7000, rprop.1, dprop.1, constrained.targets[[4]])
    res.ls <- list(res.1$X, res.2$X, res.3$X, res.4$X)
    return(res.ls)
}

res <- run.simul()

mu.i.hat <- lapply(res, colMeans)

mu.hat <- rowSums(sapply(seq(4), function(i) ws[i] * mu.i.hat[[i]]))

save(mu.hat, file=paste("results/vw/vw", idx, ".rdata", sep=""))

