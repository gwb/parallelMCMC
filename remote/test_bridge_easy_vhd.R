set.seed(123,"L'Ecuyer")

source("../bridge-sampling.R", chdir=T)
source("parallel_adaptive-hd.R")
source("../mh.R", chdir=T)
require(parallel)
require(yaml)
require(methods) # this is loaded automatically when using the R console, but not by rscript, and is necessary for the mvtnorm packages

# # # # # # # # # # # # # # # # # # # # # # # # #
# Command line parameters & Yaml configuration  #
# # # # # # # # # # # # # # # # # # # # # # # # #

# Command line
args <- commandArgs(trailingOnly=TRUE)
model.filepath <- args[1]
data.filepath <- args[2]
result.filepath <- args[3]

print(model.filepath)
print(data.filepath)

if(any(is.na(c(model.filepath, data.filepath, result.filepath)))){
    stop("One or more arguments are missing")
}

source(model.filepath)
load(data.filepath)

# yaml
params = yaml.load_file("params.yml")
sampler.ndraws <- params$bridge$sampler_ndraws
ref.ndraws <- params$bridge$ref_ndraws


get.r.centered.proposal <- function(u, sigmat){
  #fn <- function(n){rmvnorm(n, u, 2*diag(5))}
  #fn <- function(n){rmvt(n, delta=u, sigma=diag(length(u)), df=4, type="shifted")}
    fn <- function(n){rmvt(n, delta=u, sigma=sigmat, df=10, type="shifted")}
    return(fn)
}

get.d.centered.proposal <- function(u,sigmat){
    #fn <- function(x){dmvnorm(x, u, 2*diag(5))}
    #fn <- function(x){dmvt(x, delta=u, df=4, log=F)}
    fn <- function(x){dmvt(x, delta=u, sigma=sigmat, df=10, log=F)}
    return(fn)
}

K <- nrow(algo.res$centers[[1]])
ncores <- K
indx <- length(algo.res[[1]]) - 1


constrained.targets <- get.constrained.targets(K, algo.res$indicators[[indx]], dtarget)
mh.sampler.i <- get.mv.mh(sampler.ndraws, draw.normal.proposal, eval.normal.proposal, burnin=as.integer(sampler.ndraws*0.1))

res.bs <- NULL
Xs <- vector("list", K)

compute.weight.i <- function(i){
  xi <- mh.sampler.i(constrained.targets[[i]], algo.res$centers[[indx]][i,])
  u <- optim(algo.res$centers[[indx]][i,], constrained.targets[[i]], control=list(fnscale=-1))$par
  sigmat <- cov(xi)
  rq2 <- get.r.centered.proposal(u, sigmat)
  dq2 <- get.d.centered.proposal(u, sigmat)
  x2 <- rq2(ref.ndraws)
  res.i <- bridge.sampling.very.fast(constrained.targets[[i]], dq2, xi, x2)
  return(list(xi,res.i))
}

res.cwi <- mclapply(seq(K), compute.weight.i, mc.cores=ncores)
print(str(res.cwi))
#res.cwi <- lapply(seq(K), compute.weight.i)

res.bs <- do.call('c', lapply(res.cwi, function(items) items[[2]]))
Xs <- lapply(res.cwi, function(items) items[[1]])


emp.mean.ls <- lapply(Xs, colMeans)
emp.mean.mat <- do.call('rbind', emp.mean.ls)
norm.bs <- res.bs/sum(res.bs)

emp.mean <- colSums(emp.mean.mat*norm.bs)


true.mean <- colMeans(rtarget(500))

save(norm.bs, emp.mean, true.mean, file=result.filepath)
#save(norm.bs, emp.mean, true.mean, file="results/bridge_sampling_easy_vhd.rdata")
