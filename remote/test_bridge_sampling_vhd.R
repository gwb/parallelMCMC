set.seed(123,"L'Ecuyer")
#load("results/vanderwerken-1-5D-p-simple.rdata")
load("results/vanderwerken-1-5D-p-easy.rdata")
#load("results/vanderwerken-1-5D-p-easy-alt.rdata")
source("../models/vanderWerken-vhd.R")
#source("../models/vanderWerken-easy-vhd.R")
source("../bridge-sampling.R", chdir=T)
source("parallel_adaptive-hd.R")
source("../mh.R", chdir=T)

require(parallel)

get.r.centered.proposal <- function(u){
  fn <- function(n){rmvnorm(n, u, 2*diag(5))}
  return(fn)
}

get.d.centered.proposal <- function(u){
  fn <- function(x){dmvnorm(x, u, 2*diag(5))}
  return(fn)
}

K <- nrow(algo.res$centers[[1]])
ncores <- K
indx <- length(algo.res[[1]]) - 1


constrained.targets <- get.constrained.targets(K, algo.res$indicators[[indx]], dtarget)
mh.sampler.i <- get.mv.mh(4000, draw.normal.proposal, eval.normal.proposal, burnin=200)
#mh.sampler.i <- get.mv.mh(1000, draw.normal.proposal, eval.normal.proposal, burnin=200)

res.bs <- NULL
Xs <- vector("list", K)

compute.weight.i <- function(i){
  xi <- mh.sampler.i(constrained.targets[[i]], algo.res$centers[[indx]][i,])
  rq2 <- get.r.centered.proposal(algo.res$centers[[indx]][i,])
  dq2 <- get.d.centered.proposal(algo.res$centers[[indx]][i,])
  x2 <- rq2(4000)
  #x2 <- rq2(1000)
  res.i <- bridge.sampling.very.fast(constrained.targets[[i]], dq2, xi, x2)
  return(list(xi,res.i))
}

res.cwi <- mclapply(seq(K), compute.weight.i, mc.cores=ncores)
#res.cwi <- lapply(seq(K), compute.weight.i)

res.bs <- do.call('c', lapply(res.cwi, function(items) items[[2]]))
Xs <- lapply(res.cwi, function(items) items[[1]])


emp.mean.ls <- lapply(Xs, colMeans)
emp.mean.mat <- do.call('rbind', emp.mean.ls)
norm.bs <- res.bs/sum(res.bs)

emp.mean <- colSums(emp.mean.mat*norm.bs)


true.mean <- colMeans(rtarget(500))

#save(norm.bs, emp.mean, true.mean, file="results/bridge_sampling_vhd.rdata")
save(norm.bs, emp.mean, true.mean, file="results/bridge_sampling_easy_vhd.rdata")
