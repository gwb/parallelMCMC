set.seed(123,"L'Ecuyer")
#load("results/vanderwerken-1-5D-p.rdata")
load("results/vanderwerken-1-5D-p-simple.rdata")
source("models/vanderWerken-vhd.R")
source("bridge-sampling.R")
source("remote/parallel_adaptive-hd.R",chdir=T)
source("mh.R")


get.r.centered.proposal <- function(u){
  fn <- function(n){rmvnorm(n, u, diag(5))}
  return(fn)
}

get.d.centered.proposal <- function(u){
  fn <- function(x){dmvnorm(x, u, diag(5))}
  return(fn)
}

K <- 4
indx <- 6

constrained.targets <- get.constrained.targets(K, algo.res$indicators[[indx]], dtarget)
mh.sampler.i <- get.mv.mh(500, draw.normal.proposal, eval.normal.proposal, burnin=50)

res.bs <- NULL
Xs <- vector("list", K)


for(i in seq(K)){
  print(paste("Iteration", i))
  xi <- mh.sampler.i(constrained.targets[[i]], algo.res$centers[[indx]][i,])
  Xs[[i]] <- xi
  rq2 <- get.r.centered.proposal(algo.res$centers[[indx]][i,])
  dq2 <- get.d.centered.proposal(algo.res$centers[[indx]][i,])
  x2 <- rq2(500)
  #alpha.i <- geometric.bridge(constrained.targets[[i]], dq2)
  #res.i <- bridge.sampling(constrained.targets[[i]], dq2, xi, x2, alpha.i)
  res.i <- bridge.sampling.very.fast(constrained.targets[[i]], dq2, xi, x2)
  res.bs <- c(res.bs, res.i)
}


emp.mean.ls <- lapply(Xs, colMeans)
emp.mean.mat <- do.call('rbind', emp.mean.ls)
norm.bs <- res.bs/sum(res.bs)

emp.mean <- colSums(emp.mean.mat*norm.bs)


true.mean <- colMeans(rtarget(500))
