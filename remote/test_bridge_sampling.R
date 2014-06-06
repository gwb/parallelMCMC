
# The purpose here is to test different implementations of bridge sampling, to speed things
# up. This results in a "very fast" implementation that divides the execution time by a little
# of 7, compared to the original version. The bottleneck is in the evaluation of mvnorm. We can't
# do away with the "apply" since the first argument, which will typically be the constrained target
# is not vectorized (and vectorization here would be superficial and won't speed things up)

# Obviously, I tested that all the implementations give exactly the same results.


source("bridge-sampling.R")
require("mvtnorm")

rq1 <- function(n) {rmvnorm(n,c(1,2), matrix(c(1,0,0,1), byrow=T, nrow=2))}
dq1 <- function(x) {dmvnorm(x,c(1,2), matrix(c(1,0,0,1), byrow=T, nrow=2))}
rq2 <- function(n) {rmvnorm(n,c(3,2), matrix(c(1,0,0,1), byrow=T, nrow=2))}
dq2 <- function(x) {dmvnorm(x,c(3,2), matrix(c(1,0,0,1), byrow=T, nrow=2))}

alpha <- geometric.bridge(dq1, dq2)

res.is <- NULL
for(i in seq(15)){
  tmp.res <- importance.sampling(dq1, rq2, dq2, 1000)
  res.is <- c(res.is, tmp.res)
}
mean(res.is)
var(res.is)

set.seed(123)
foo1 <- function(){
alpha <- geometric.bridge(dq1, dq2)
res.bs <- NULL
for(i in seq(4)){
  x1 <- rq1(500)
  x2 <- rq2(500)
  tmp.res <- bridge.sampling(dq1, dq2, x1, x2, alpha)
  res.bs <- c(res.bs, tmp.res)
}
return(c(mean(res.bs),var(res.bs)))
}

set.seed(123)
foo2 <- function(){
alpha.fast <- geometric.bridge.fast(dq1, dq2)
res.bs <- NULL
for(i in seq(4)){
  x1 <- rq1(500)
  x2 <- rq2(500)
  tmp.res <- bridge.sampling.fast(dq1, dq2, x1, x2, alpha.fast)
  res.bs <- c(res.bs, tmp.res)
}
return(c(mean(res.bs),var(res.bs)))
}

set.seed(123)
foo3 <- function(){
res.bs <- NULL
for(i in seq(4)){
  x1 <- rq1(500)
  x2 <- rq2(500)
  tmp.res <- bridge.sampling.very.fast(dq1, dq2, x1, x2)
  res.bs <- c(res.bs, tmp.res)
}
return(c(mean(res.bs),var(res.bs)))
}

set.seed(123)
Rprof("bridge-slow.rprof")
replicate(10,foo1())
Rprof(NULL)
rslow <- summaryRprof("bridge-slow.rprof")

set.seed(123)
Rprof("bridge-fast.rprof")
replicate(10,foo2())
Rprof(NULL)
rfast <- summaryRprof("bridge-fast.rprof")

set.seed(123)
Rprof("bridge-very-fast.rprof")
replicate(10,foo3())
Rprof(NULL)
rveryfast <- summaryRprof("bridge-very-fast.rprof")

