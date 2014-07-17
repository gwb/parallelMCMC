set.seed(123,"L'Ecuyer")

source("../graphics.R")
source("../models/vanderWerken-easy-vhd.R")
load("results/vanderwerken-1-5D-p-easy.rdata")


# reference
ref <- rtarget.clust(200)

ref.draws <- ref[[1]]
ref.clust <- ref[[2]]

bla <- plot.vhd.clusters(ref.draws, clust=ref.clust)


# comparing partitioning at last iteration with reference partitionin

## extracting the cluster indicator function found by our algorithm
indx <- length(algo.res$indicators)
clust.fn <- algo.res$indicators[[indx]]

## preprocessing
algo.dt <- plot.vhd.clusters(ref.draws, cluster.fn=clust.fn)
dt <- rbind(algo.dt[[1]], algo.dt[[1]])
dt$clust <- c(rep(ref.clust, 5), algo.dt[[1]]$clust)
dt$method <- rep(c("true", "spectral"), each=length(algo.dt[[1]][,1]))

## plotting
g <- ggplot(data=dt, aes(x,y)) + geom_point(aes(color=factor(clust))) + facet_grid(method~proj)
ggsave(filename="results/algo_vs_ref_partitioning.pdf", plot=g, width=16, height=10)


# plotting the progression of the partitioning over multiple iterations of the algo

algo.dt <- plot.multiple.vhd.clusters(ref.draws, cluster.fn.ls = algo.res$indicators, centers.ls = algo.res$centers)

g <- ggplot(data=algo.dt[[1]], aes(x,y)) + geom_point(aes(color=factor(clust))) + facet_grid(ind~proj) +
  geom_point(data=algo.dt[[2]], aes(x,y))

ggsave(filename="results/partitioning_progression.pdf", plot=g, width=16, height=10)
