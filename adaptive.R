source("mh.R")
source("approx.transitions.R")
source("spectral-methods.R")
source("graphics.R")


set.seed(2)

# # # # # # 
#  SETUP  #
# # # # # #

K <- 2 # number of components of mixture
s <- 0.15 # standard dev of mixture components
inter <- 2 # interval between means 

prop.sd <- 0.15 # standard dev of normal proposal
mcmc.iter <- 5000

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


# # # # # # # # # # # # # # # #
# RUNNING INDEPENDENT CHAINS  #
# # # # # # # # # # # # # # # #

all.res <- NULL
for( i in seq(K) ){  
  mh.res <- run.mh(inter*i, mcmc.iter, draw.normal.proposal, eval.normal.proposal, dtarget) 
  #all.draws <- c(all.draws, mh.res[[1]][,1])
  all.res <- rbind(all.res, mh.res[[1]])
}

X <- all.res[,1] # all the states
L <- all.res[,2] # all the proposals
a <- ifelse(all.res[,3] <= 1, all.res[,3], 1) # the acceptance function

# # # # # # # # # # # # # # # # # # # #
# APPROXIMATING THE TRANSITION KERNEL #
# # # # # # # # # # # # # # # # # # # #

# updates min and max vals to make sure that they capture all
# values of X and L
min.val <- min(c(min.val, X, L)) - 0.05
max.val <- max(c(max.val, X, L)) + 0.05

discretization <- seq(min.val, max.val, length.out=num.bins)
intervals <- levels(cut(seq(min.val,max.val), discretization))

K.approx <- get.K.approx(discretization, X, L, a)

clusters <- asym.spectral.clustering(K.approx, 2)




# # # # # # # # # 
# SOME PLOTTING #
# # # # # # # # #

pdf(file="normal-mixture-k-2.pdf")
curve(dtarget, from=min.val, to = max.val, n=5000, add=FALSE, col="blue", ylab=NA)
hist(X, breaks=1000, col="red", freq=FALSE, add=TRUE)
curve(dtarget, from=min.val, to = max.val, n=5000, add=TRUE, col="blue", lwd=2)
clusters.intervals <- cbind(get.numerical.intervals(intervals), clusters)
my.cols <- c(rgb(1,0,0,alpha=0.2), rgb(0,1,0, alpha=0.2))
add.rects(clusters.intervals, my.cols)
title(main="Target density \n with histogram from independent MH and calculated partition",
      ylab="density")
dev.off()


# # # # # # # # #
# SAVING STATE  #
# # # # # # # # #

save(K, s, inter, prop.sd, mcmc.iter, min.val, max.val, num.bins,
     dtarget,
     draw.normal.proposal,
     eval.normal.proposal,
     X, L, a,
     discretization, intervals, K.approx, clusters,
     file="normal-mixture-k-2.Rdata")
