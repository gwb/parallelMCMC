require(parallel)
source("../mh.R")
source("../clustering.R")


# args: init.poins  - matrix where each row is an initial point
#       mcmcsampler - a sampler that takes an initial points, and returns
#                     a matrix of draws where each row is a draw
init.exploration.test <- function(dtarget, mcmc.sampler, init.points){
  if(is.vector(init.points)){
    init.points <- as.matrix(init.points)
  }

  #init.points <- init.points[c(1,2),]
  
  res <- NULL
  for(i in seq(nrow(init.points))){
    draws <- mcmc.sampler(dtarget, init.points[i,])
    res <- rbind(res, draws)
  }

  return(res)
}

# args: init.poins  - matrix where each row is an initial point
#       mcmcsampler - a sampler that takes an initial points, and returns
#                     a matrix of draws where each row is a draw
init.exploration.old <- function(dtarget, mcmc.sampler, init.points){
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



init.exploration <- function(dtarget, mcmc.sampler, init.points){
  if(is.vector(init.points)){
    init.points <- as.matrix(init.points)
  }

  draws <- do.call('rbind', mclapply(seq(nrow(init.points)), function(i) mcmc.sampler(dtarget, init.points[i,])))

  return(draws)
}


