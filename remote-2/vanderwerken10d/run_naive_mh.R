require(methods)
require(mvtnorm)
source("../../mh.R", chdir=T)
source("../../clustering.R", chdir=T)
source("../../bridge-sampling.R", chdir=T)

args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[1])

load('results/exploration.rdata')

res.1 <- run.mv.mh(u.ls[[1]], 10000, rprop.1, dprop.1, dtarget)
res.2 <- run.mv.mh(u.ls[[2]], 10000, rprop.1, dprop.1, dtarget)
res.3 <- run.mv.mh(u.ls[[3]], 10000, rprop.1, dprop.1, dtarget)
res.4 <- run.mv.mh(u.ls[[4]], 10000, rprop.1, dprop.1, dtarget)

res <- rbind(res.1$X, res.2$X, res.3$X, res.4$X)

mu.hat <- colMeans(res)


save(mu.hat, file=paste("results/naive/naive", idx, ".rdata", sep=""))
