require(methods)
require(mvtnorm)
require(parallel)
source("../../mh.R", chdir=T)
source("../../clustering.R", chdir=T)

args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[1])

load("results/exploration-vc.rdata")

get.clust.init <- function(X){
    clust.1.idx <- which(vc$cluster==1)
    clust.2.idx <- which(vc$cluster==2)
    return(list(X[sample(clust.1.idx,1),], X[sample(clust.2.idx,1),]))
}

run.simul <- function(centers){
    res.vc.1 <- run.mv.mh(vc.centers[[1]], 15000, rprop.1, dprop.1, vc.constrained.targets[[1]])
    res.vc.2 <- run.mv.mh(vc.centers[[2]], 15000, rprop.1, dprop.1, vc.constrained.targets[[2]])
    return(ws.vc[1] * colMeans(res.vc.1$X) + ws.vc[2] * colMeans(res.vc.2$X))
}

centers <- get.clust.init(xs)

res <- run.simul(centers)

#res <- do.call('cbind', mclapply(seq(64), run.simul, mc.cores = 32))

save(res, centers, file=paste("results/vc/vc", idx, ".rdata", sep=""))
#save(res, centers, file="results/simul-vc.rdata")
