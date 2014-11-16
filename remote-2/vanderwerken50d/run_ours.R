require(methods)
require(mvtnorm)
#source("../../parallel-tempering.R", chdir=T)
source("../../mh.R", chdir=T)
#source("../../clustering.R", chdir=T)
#source("../../bridge-sampling.R", chdir=T)

args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[1])

load('results/exploration.rdata')

true.ws <- c(0.1,0.2,0.3,0.4)
true.mean <- rowSums(do.call('cbind',lapply(seq(4), function(i) true.ws[i] * u.ls[[i]])))


c.run.simul <- function(){
    res.1 <- run.mv.mh.fast(vc$centers[1,], 10000, rprop, dprop, constrained.targets[[1]])
    res.2 <- run.mv.mh.fast(vc$centers[2,], 10000, rprop, dprop, constrained.targets[[2]])
    res.3 <- run.mv.mh.fast(vc$centers[3,], 10000, rprop, dprop, constrained.targets[[3]])
    res.4 <- run.mv.mh.fast(vc$centers[4,], 10000, rprop, dprop, constrained.targets[[4]])
    res.ls <- list(res.1$X, res.2$X, res.3$X, res.4$X)
    return(res.ls)
}

res <- c.run.simul()

mu.i.hat <- lapply(res, colMeans)

mu.hat <- rowSums(sapply(seq(4), function(i) ws[i] * mu.i.hat[[i]]))

save(mu.hat, file=paste("results/ours/ours", idx, ".rdata", sep=""))
#
#
#
#mu.hat.vc <- rowSums(sapply(seq(4), function(i) ws[i] * vc$centers[i,]))
#
#
#alt.true.mean <- ws[3]*u.ls[[1]] + ws[1]*u.ls[[2]] + ws[2]*u.ls[[3]] + ws[4]*u.ls[[4]]
#
#
#
#sum((mu.i.hat[[3]] - u.ls[[1]])^2)
#sum((vc$centers[3,] - u.ls[[1]])^2)
#
#
#sum((mu.i.hat[[1]] - u.ls[[2]])^2)
#sum((vc$centers[1,] - u.ls[[2]])^2)
#
#sum((mu.i.hat[[2]] - u.ls[[3]])^2)
#sum((vc$centers[2,] - u.ls[[3]])^2)
#
#sum((mu.i.hat[[4]] - u.ls[[4]])^2)
#sum((vc$centers[4,] - u.ls[[4]])^2)
#
#
#
#
