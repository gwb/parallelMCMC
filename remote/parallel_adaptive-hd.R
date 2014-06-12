require(parallel)
source("../mh.R")
source("../clustering.R")


init.exploration <- function(dtarget, mcmc.sampler, init.points,ncores=2){
  if(is.vector(init.points)){
    init.points <- as.matrix(init.points)
  }

  draws <- do.call('rbind', mclapply(seq(nrow(init.points)),
                                     function(i) mcmc.sampler(dtarget, init.points[i,]), mc.cores=ncores))

  save(draws, file="output_p.rdata")
  
  return(draws)
}


# TODO: put parallel back
adaptive.sampling <- function(dtarget, mcmc.sampler, n.clust, init.draws, clust.method, algo.iter, ncores=2){

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
  #res.draws <- vector('list', algo.iter+1)
  #res.draws[[1]] <- init.draws
  
  for(i in seq(algo.iter)){
    print(paste("Iteration:", i))
    all.draws <- NULL
    
    # parallel mcmc
    #for(j in seq(n.clust)){
    #  print(paste("... Chain:", j))
    #  draws <- mcmc.sampler(constrained.targets[[j]], clust.centers[j,])
    #  all.draws <- rbind(all.draws, draws)
    #}

    # parallel mcmc
    all.draws <- do.call('rbind',
                         mclapply(seq(n.clust),
                                  function(j) mcmc.sampler(constrained.targets[[j]], clust.centers[j,]),
                                  mc.cores=ncores))
    
    # clustering
    clust.res <- clust.method(all.draws, n.clust) 
    clust.indicator <- clust.res$indicator
    clust.centers <- clust.res$centers

    # transform the target
    constrained.targets <- get.constrained.targets(n.clust, clust.indicator, dtarget)

    # save cluster indicator
    clust.indicator.list[[i+1]] <- clust.indicator
    clust.centers.list[[i+1]] <- clust.centers

    # save draws
    #res.draws[[i+1]] <- all.draws
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
