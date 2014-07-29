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
params = yaml.load_file("params.yml")
explorer.ndraws <- params$adaptive$explorer_ndraws
sampler.ndraws <- params$adaptive$sampler_ndraws
nsub <- params$adaptive$nsub
niter <- params$adaptive$niter


# # # # # # # # # # # #
# Running the sampler
# # # # # # # # # # # #

# We initialize at the actual centers. See how we do in this simple case
#init.points <- do.call('rbind', lapply(1:4, function(i) do.call('rbind', u.ls)))
init.points <- do.call('rbind', u.ls)
save(init.points, file="debug-exploration.rdata")

# Number of desired partitions
K <- 4

# create the samplers
mh.sampler <- get.mv.mh(sampler.ndraws, draw.normal.proposal, eval.normal.proposal, burnin=as.integer(sampler.ndraws*0.1))
mh.explorer <- get.mv.mh(explorer.ndraws, draw.normal.proposal, eval.normal.proposal, burnin=as.integer(explorer.ndraws*0.1))

# explore the surface
print("BEGIN: exploration")
init.draws <- init.exploration(dtarget, mh.explorer, init.points, ncores=ncores)
print("END: exploration")

# setting up clustering method
spec.clust <- function(X, C){
    #return(spectral.clustering(X, C, L, nsub=100,L.vec=L.vec))
    return(spectral.clustering(X, C, L.vec.sym, get.K.hat, .get.closeness.centers, nsub=nsub))
}

# runs the adaptive sampling algorithm
algo.res <- adaptive.sampling(dtarget, mh.sampler, K, init.draws, spec.clust, niter, ncores=K)



    
#save(u.ls, sig.ls, init.draws, init.points, algo.res, file="results/vanderwerken-1-5D-p-easy.rdata")
save(u.ls, sig.ls, init.draws, init.points, algo.res, model.filepath, params, file=output.filepath)

