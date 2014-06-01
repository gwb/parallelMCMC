source("mh.R")
source("clustering.R")

# args: init.poins  - matrix where each row is an initial point
#       mcmcsampler - a sampler that takes an initial points, and returns
#                     a matrix of draws where each row is a draw
init.exploration <- function(dtarget, mcmc.sampler, init.points){
  if(is.vector(init.points)){
    init.points <- as.matrix(init.points)
  }

  res <- NULL
  for(i in seq(nrow(init.points))){
    draws <- mcmc.sampler(dtarget, init.points[i,])
    res <- rbind(res, draws)
  }

  return(res)
}


adaptive.sampling <- function(dtarget, mcmc.sampler, n.clust, init.draws, clust.method, algo.iter){

  # start by obtaining initial clusters (from initialization data)
  clust.res <- clust.method(init.draws, n.clust)
  clust.indicator <- clust.res$indicator
  clust.centers <- clust.res$centers # matrix where each row is a the center of a cluster

  # transforms the target
  constrained.targets <- get.constrained.targets(n.clust, clust.indicator, dtarget)


  # stuff we want to return at the end
  clust.indicator.list <- vector('list', algo.iter+1)
  clust.indicator.list[[1]] <- clust.indicator
  clust.centers.list <- vector('list', algo.iter+1)
  clust.centers.list[[1]] <- clust.centers

  
  for(i in seq(algo.iter)){
    print(paste("Iteration:", i))
    all.draws <- NULL
    
    # parallel mcmc
    for(j in seq(n.clust)){
      print(paste("... Chain:", j))
      draws <- mcmc.sampler(constrained.targets[[j]], clust.centers[j,])
      all.draws <- rbind(all.draws, draws)
    }

    # clustering
    clust.res <- clust.method(all.draws, n.clust) 
    clust.indicator <- clust.res$indicator
    clust.centers <- clust.res$centers

    # transform the target
    constrained.targets <- get.constrained.targets(n.clust, clust.indicator, dtarget)

    # save cluster indicator
    clust.indicator.list[[i+1]] <- clust.indicator
    clust.centers.list[[i+1]] <- clust.centers
  }

  return(list(indicators = clust.indicator.list, centers = clust.centers.list))

}


# # # # 
# Different samplers
# # # #

get.mv.mh <- function(niter, rproposal, dproposal, burnin=0){
  fn <- function(dtarget, x0){
    force(x0)
    return(run.mv.mh(x0, niter, rproposal, dproposal, dtarget, burnin)$X)
  }
  return(fn)
}
