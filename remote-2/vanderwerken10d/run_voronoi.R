require(methods)
require(mvtnorm)
source("../../parallel-tempering.R", chdir=T)
source("../../mh.R", chdir=T)
source("../../clustering.R", chdir=T)
source("../../bridge-sampling.R", chdir=T)

args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[1])

load('results/exploration-vc.rdata')

true.mean <- 0.1 * u.ls[[1]] + 0.2 * u.ls[[2]] + 0.3 * u.ls[[3]] + 0.4 * u.ls[[4]]

get.clust.init <- function(X){
    clust.1.idx <- which(vc$cluster==1)
    clust.2.idx <- which(vc$cluster==2)
    clust.3.idx <- which(vc$cluster==3)
    clust.4.idx <- which(vc$cluster==4)
    return(list(X[sample(clust.1.idx,1),], X[sample(clust.2.idx,1),],
                X[sample(clust.3.idx,1),], X[sample(clust.4.idx,1),]))
}



run.simul <- function(centers){
    res.1 <- run.mv.mh(centers[[1]], 5000, rprop.1, dprop.1, vc.constrained.targets[[1]])
    res.2 <- run.mv.mh(centers[[2]], 5000, rprop.1, dprop.1, vc.constrained.targets[[2]])
    res.3 <- run.mv.mh(centers[[3]], 5000, rprop.1, dprop.1, vc.constrained.targets[[3]])
    res.4 <- run.mv.mh(centers[[4]], 5000, rprop.1, dprop.1, vc.constrained.targets[[4]])
    return(ws.vc[1] * colMeans(res.1$X) +
           ws.vc[2] * colMeans(res.2$X) +
           ws.vc[3] * colMeans(res.3$X) +
           ws.vc[4] * colMeans(res.4$X))
}

centers <- get.clust.init(xs)

mu.hat <- run.simul(centers)

#mu.i.hat <- lapply(res, colMeans)

#mu.hat <- rowSums(sapply(seq(4), function(i) ws[i] * mu.i.hat[[i]]))

save(mu.hat, file=paste("results/vc/vc", idx, ".rdata", sep=""))

