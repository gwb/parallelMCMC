source("mh.R")
source("approx.transitions.R")
source("spectral-methods.R")
source("graphics.R")


set.seed(2)


# # # # # # 
#  SETUP  #
# # # # # #

M <- 2 # number of // chains
K <- 2 # number of components of mixture
C <- 2 # number of clusters to look for
s <- 0.15 # standard dev of mixture components
inter <- 2 # interval between means 

prop.sd <- 0.1 # standard dev of normal proposal  # TRY WITH BADLY CALIBRATED
mcmc.iter <- 500 # very small number of iterations
x0 <- c(1,5) # initial values of samplers

min.val <- 1 # minimum value of interest for plotting and discretization
max.val <- 5 # maximum value of interest for ..
num.bins <- 101 # number of bins for discretization


# defining the multimodal target 
dtarget <- function(x){
  return((1/K) * sum(sapply(seq(1,K), function(i) dnorm(x, mean=inter*i, sd=s))))
}
dtarget <- Vectorize(dtarget)

# function to draw from the proposal distribution
draw.normal.proposal <- function(x){
  return(rnorm(1, mean=x, sd=prop.sd))
}

# function to evaluate the transition between old state and new state
eval.normal.proposal <- function(new.state, old.state){
  return(dnorm(new.state, mean=old.state, sd=prop.sd))
}


# Simulates the adaptive sampling algorithm in 1D, using spectral clustering,
# and plots the partition found at teach iteration of the algorithm (plotting
# is done to a file, in pdf)
# Args - x0        : a vector
#        n.rounds  : number of rounds of the algorithm to run
#        M         : number of parallel chains to run
#        C         : number of clusters to look for
#        mcmc.iter : number of iterations of the mcmc algorithm
run.adaptive <- function(x0, n.rounds, M=M, C=C, mcmc.iter=mcmc.iter){

  for(j in seq(n.rounds)){
    print(paste("Current round: ", j, sep=""))
    cat("Initial values for chains: ", x0, "\n")
    all.res <- NULL
    for(i in seq(M)){
      mh.res <- run.mh(x0[i], mcmc.iter, draw.normal.proposal, eval.normal.proposal, dtarget) 
      all.res <- rbind(all.res, mh.res[[1]])
    }

    X <- all.res[,1] # all the states
    L <- all.res[,2] # all the proposals
    a <- ifelse(all.res[,3] <= 1, all.res[,3],1) # the acceptance function
    
    c.min.val <- min(c(min.val, X, L)) - 0.05
    c.max.val <- max(c(max.val, X, L)) + 0.05
    discretization <- seq(c.min.val, c.max.val, length.out=num.bins)
    intervals <- levels(cut(seq(c.min.val,c.max.val), discretization))
    #K.approx.res <- get.K.approx.new.plus(discretization, X, L, a)
    K.approx.res <- get.K.approx.new(discretization, X, L, a)
    K.approx <- K.approx.res[[1]]

    clusters <- asym.spectral.clustering(K.approx, C)
    clusters.post <- post.process.clusters(clusters, K.approx.res[[2]], C)
    K.partition.info <- K.approx.res[[2]]
    save(K.approx, clusters, clusters.post, K.partition.info, file=paste("K_and_clusters_",j,".rdata",sep=""))
    
    clusters <- clusters.post
    
    # plotting things
    {
      pdf(paste("partition_round_",j,".pdf", sep=""))
      curve(dtarget, from=c.min.val, to = c.max.val, n=5000, add=FALSE, col="blue", ylab=NA)
      hist(X, breaks=1000, col="red", freq=FALSE, add=TRUE)
      curve(dtarget, from=c.min.val, to = c.max.val, n=5000, add=TRUE, col="blue", lwd=2)
      clusters.intervals <- cbind(get.numerical.intervals(intervals), clusters)
      add.rects(clusters.intervals, my.cols)
      for(x in x0){
        abline(v=x, lwd=2, col="yellow")
      }
      dev.off()
    }

    # updating center of clusters
    x0 = rep(0, C)
    for(k in seq(C)){
      x0[k] <- discretization[as.integer(median(which(clusters==k)))]
    }
    
  } 
}


# This function differs from run.adaptive in that it uses not just the draws from the current
# round, but the draws from all previous rounds too, to compute the new partition. Memory refers
# to this use of past information. The parameters are otherwise the same
run.adaptive.memory <- function(x0, n.rounds, M=M, K=K, C=C, mcmc.iter=mcmc.iter){

  all.res <- NULL
  for(j in seq(n.rounds)){
    print(paste("Current round: ", j, sep=""))
    cat("Initial values for chains: ", x0, "\n")
    for(i in seq(M)){
      mh.res <- run.mh(x0[i], mcmc.iter, draw.normal.proposal, eval.normal.proposal, dtarget) 
      all.res <- rbind(all.res, mh.res[[1]])
    }

    X <- all.res[,1] # all the states
    L <- all.res[,2] # all the proposals
    a <- ifelse(all.res[,3] <= 1, all.res[,3],1) # the acceptance function
    
    c.min.val <- min(c(min.val, X, L)) - 0.05
    c.max.val <- max(c(max.val, X, L)) + 0.05
    discretization <- seq(c.min.val, c.max.val, length.out=num.bins)
    intervals <- levels(cut(seq(c.min.val,c.max.val), discretization))
    K.approx.res <- get.K.approx.new.plus(discretization, X, L, a)
    #K.approx.res <- get.K.approx.new(discretization, X, L, a)
    K.approx <- K.approx.res[[1]]

    clusters <- asym.spectral.clustering(K.approx, C)
    clusters <- post.process.clusters(clusters, K.approx.res[[2]], C)

    # plotting things
    {
      pdf(paste("partition_round_",j,".pdf", sep=""))
      curve(dtarget, from=c.min.val, to = c.max.val, n=5000, add=FALSE, col="blue", ylab=NA)
      hist(X, breaks=1000, col="red", freq=FALSE, add=TRUE)
      curve(dtarget, from=c.min.val, to = c.max.val, n=5000, add=TRUE, col="blue", lwd=2)
      clusters.intervals <- cbind(get.numerical.intervals(intervals), clusters)
      add.rects(clusters.intervals, my.cols)
      for(x in x0){
        abline(v=x, lwd=2, col="yellow")
      }
      dev.off()
    }

    # updating center of clusters
    x0 = rep(0, C)
    for(k in seq(C)){
      x0[k] <- discretization[as.integer(median(which(clusters==k)))]
    }
    
  } 
}


run.adaptive(x0, 5, M=M, K=K, C=C, mcmc.iter=mcmc.iter)



M <- 3 # number of // chains
K <- 3 # number of components of mixture
C <- 3 # number of clusters to look for
s <- 0.1 # standard dev of mixture components
inter <- 3 # interval between means 

prop.sd <- 0.01 # standard dev of normal proposal  # TRY WITH BADLY CALIBRATED
mcmc.iter <- 2000 # very small number of iterations
x0 <- c(1,3,5) # initial values of samplers

min.val <- 1 # minimum value of interest for plotting and discretization
max.val <- 12 # maximum value of interest for ..
num.bins <- 50 # number of bins for discretization


# defining the multimodal target 
dtarget <- function(x){
  return((1/K) * sum(sapply(seq(1,K), function(i) dnorm(x, mean=inter*i, sd=s))))
}
dtarget <- Vectorize(dtarget)

# function to draw from the proposal distribution
draw.normal.proposal <- function(x){
  return(rnorm(1, mean=x, sd=prop.sd))
}

# function to evaluate the transition between old state and new state
eval.normal.proposal <- function(new.state, old.state){
  return(dnorm(new.state, mean=old.state, sd=prop.sd))
}


curve(dtarget, from=min.val, to = max.val, n=5000)

run.adaptive(x0, 10, M=M, K=K, C=C, mcmc.iter=mcmc.iter)

run.adaptive.memory(x0, 30, M=M, K=K, C=C, mcmc.iter=mcmc.iter)


load("K_and_clusters_8.rdata")

