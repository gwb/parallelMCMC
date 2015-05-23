source("model-univariate.R")
source("performance.R") 
source("../mh.R")

args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[1])

# list has length 39
mu.ls <- seq(0.2, 4, 0.1)

mu <- mu.ls[idx]
tau <- 0.5
dproposal <- get.dproposal(tau)
rproposal <- get.rproposal(tau)
dtarget <- get.dtarget(mu)
dtransition <- get.dtransition(tau, mu, dproposal, dtarget)

nsimul <- 10
midpoint.res.ls <- vector('list', nsimul)
sg.list <- vector('list', nsimul)
for(i in seq(nsimul)){
    midpoint.res <- get.midpoint(6000, 400, 1000, rproposal, dproposal, dtarget)
    Yi <- as.vector(midpoint.res[[1]])
    midpoint <- midpoint.res[[2]]
    clust <- midpoint.res[[4]][[1]]
    
    empirical.sg <- min(get.empirical.spectral.gaps(Yi, midpoint, dtransition))
    exact.sg <- min(get.exact.spectral.gaps(mu, midpoint, dtransition, 300))

    midpoint.res.ls[[i]] <- midpoint.res
    sg.list[[i]] <- c(empirical.sg, exact.sg)
}

save(mu, midpoint.res.ls, sg.list,
     file=paste("results/res_", idx, sep=""))
