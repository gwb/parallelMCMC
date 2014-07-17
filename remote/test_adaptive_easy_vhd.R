source("parallel_adaptive-hd.R")
source("../bridge-sampling.R")
source("../graphics.R")
source("../clustering.R")
require("yaml")

set.seed(123, "L'Ecuyer")

# # # # # # # # # # # # # # # # # # # # # # # # #
# Command line parameters & Yaml configuration  #
# # # # # # # # # # # # # # # # # # # # # # # # #

# Command line
args <- commandArgs(trailingOnly=TRUE)
ncores <- as.numeric(args[1])
model.filepath <- args[2] 
output.filepath <- args[3]


if(is.na(model.filepath)){
    model.filepath <- "../models/vanderWerken-easy-vhd.R"
}

if(is.na(output.filepath)){
    output.filepath <- "results/vanderwerken-1-5D-p-easy.rdata"
}

#source("../models/vanderWerken-easy-vhd.R")
source(model.filepath)

# Yaml


# # # # # # # # # # # #
# Running the sampler
# # # # # # # # # # # #

# We initialize at the actual centers. See how we do in this simple case
init.points <- do.call('rbind', lapply(1:4, function(i) do.call('rbind', u.ls)))


# Number of desired partitions
K <- 4

# create the samplers
mh.sampler <- get.mv.mh(1000, draw.normal.proposal, eval.normal.proposal, burnin=200)
mh.explorer <- get.mv.mh(4000, draw.normal.proposal, eval.normal.proposal, burnin=1000)

# explore the surface
init.draws <- init.exploration(dtarget, mh.explorer, init.points, ncores=ncores)

# setting up clustering method
spec.clust <- function(X, C){
    #return(spectral.clustering(X, C, L, nsub=100,L.vec=L.vec))
    return(spectral.clustering(X, C, L.vec.sym, get.K.hat, .get.closeness.centers, nsub=200))
}

# runs the adaptive sampling algorithm
algo.res <- adaptive.sampling(dtarget, mh.sampler, K, init.draws, spec.clust, 7, ncores=K)



    
#save(u.ls, sig.ls, init.draws, init.points, algo.res, file="results/vanderwerken-1-5D-p-easy.rdata")
save(u.ls, sig.ls, init.draws, init.points, algo.res, file=output.filepath)

