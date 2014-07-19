set.seed(123,"L'Ecuyer")

source("../graphics.R")


# # # # # # # # # # # # # # # # # # # # # # # # #
# Command line parameters & Yaml configuration  #
# # # # # # # # # # # # # # # # # # # # # # # # #

# Command line
args <- commandArgs(trailingOnly=TRUE)
model.filepath <- args[1] 
data.filepath <- args[2]
output.dir <- args[3]

if(is.na(model.filepath)){
    model.filepath <- "../models/vanderWerken-easy-vhd.R"
}

if(is.na(data.filepath)){
    data.filepath <- "results/vanderwerken-1-5D-p-easy.rdata"
}

if(is.na(output.dir)){
    output.dir <- "results/"
}

print(model.filepath)
print(data.filepath)
print(output.dir)


load(data.filepath)
source(model.filepath)


# # # # # # # # # # # #
# Do the actual plotting
# # # # # # # # # # # #

# reference
ref <- rtarget.clust(200)

ref.draws <- ref[[1]]
ref.clust <- ref[[2]]

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
ggsave(filename=paste(output.dir, "/algo_vs_ref_partitioning.pdf", sep=""),
       plot=g, width=16, height=10)


# plotting the progression of the partitioning over multiple iterations of the algo

algo.dt <- plot.multiple.vhd.clusters(ref.draws, cluster.fn.ls = algo.res$indicators, centers.ls = algo.res$centers)

g <- ggplot(data=algo.dt[[1]], aes(x,y)) + geom_point(aes(color=factor(clust))) + facet_grid(ind~proj) +
  geom_point(data=algo.dt[[2]], aes(x,y))

ggsave(filename=paste(output.dir, "/partitioning_progression.pdf", sep=""),
       plot=g, width=16, height=10)
