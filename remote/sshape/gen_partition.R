require(methods)
require(mvtnorm)
source("../../mh.R", chdir=T)
source("../../clustering.R", chdir=T)
source("../../bridge-sampling.R", chdir=T)

set.seed(123)

# # # # #
# Model #
# # # # #

rshape1 <- function(){
    thet <- runif(1, -4*pi/6, 4*pi/6)
    r <- runif(1,1,1.1)
    return(c(0,-1) + c(r * cos(thet), r * sin(thet)))
}

rshape2 <- function(){
    thet <- runif(1, 2*pi/6, 10*pi/6)
    r <- runif(1,1,1.1)
    return(c(0, 0) + c(r * cos(thet), r * sin(thet)))
}

rfn <- function(){
    idx <- sample(c(1,2), 1)
    rshapes <- list(rshape1, rshape2)
    return(rshapes[[idx]]())
}


isshape1 <- function(x){
    x <- as.vector(x) - c(0,-1)
    norm.x <- sqrt(sum(x^2))
    cosa <- as.vector(t(x) %*% c(1,0)) / norm.x
    if(cosa > cos(4*pi/6)){
        if(norm.x > 1 & norm.x < 1.1){
            return(TRUE)
        } else {
            return(FALSE)
        }
    }else{
        return(FALSE)
    }
}

isshape2 <- function(x){
    x <- as.vector(x)
    norm.x <- sqrt(sum(x^2))
    cosa <- as.vector(t(x) %*% c(1,0)) / norm.x
    if(cosa < cos(2*pi/6)){
        if(norm.x > 1 & norm.x < 1.1){
            return(TRUE)
        } else {
            return(FALSE)
        }
    } else {
        return(FALSE)
    }
}

dfn <- function(x){
    if(isshape1(x) || isshape2(x)){
        return(1 / ( 2 * 2/3 * pi * (1.1^2 - 1)))
    } else {
        return(0)
    }
}

rprop.1 <- function(ostate){
    return(rmvnorm(1, ostate, 0.015 * diag(2)))
}

dprop.1 <- function(nstate, ostate){
    return(dmvnorm(nstate, ostate, 0.015 * diag(2)))
}


# # # # # # # # 
# Exploration #
# # # # # # # #

res1 <- run.mv.mh(c(1.05,-1), 5000, rprop.1, dprop.1, dfn)
res2 <- run.mv.mh(c(-1.05,0), 5000, rprop.1, dprop.1, dfn)

xs <- rbind(res1$X, res2$X)

#plot(xs[,2] ~ xs[,1], xlim=c(-2,2), ylim=c(-2,2))

# # # # # # # # # # # #
# Spectral Clustering
# # # # # # # # # # # #

sc <- spectral.clustering(xs, 2, dprop.1, get.K.hat, nsub=700)
sc.indicator.fn <- sc$indicator

sc.constrained.targets <- get.constrained.targets(2, sc.indicator.fn, dfn)

res.1.sc <- run.mv.mh(sc$centers[1,], 10000, rprop.1, dprop.1, sc.constrained.targets[[1]])
res.2.sc <- run.mv.mh(sc$centers[2,], 10000, rprop.1, dprop.1, sc.constrained.targets[[2]])

# # # # # # # # # # # #
# Voronoi Clustering  #
# # # # # # # # # # # #

vc <- voronoi.clustering(xs, 2)
vc.indicator.fn <- vc$indicator

vc.constrained.targets <- get.constrained.targets(2, vc.indicator.fn, dfn)

get.clust.init <- function(X){
    clust.1.idx <- which(vc$cluster==1)
    clust.2.idx <- which(vc$cluster==2)
    return(list(X[sample(clust.1.idx,1),], X[sample(clust.2.idx,1),]))
}

vc.centers <- get.clust.init(xs)

res.1.vc <- run.mv.mh(vc.centers[[1]], 10000, rprop.1, dprop.1, vc.constrained.targets[[1]])
res.2.vc <- run.mv.mh(vc.centers[[2]], 10000, rprop.1, dprop.1, vc.constrained.targets[[2]])

# # # # # # # # # #
# Bridge Sampling #
# # # # # # # # # #


get.r.centered.proposal <- function(u, sigmat){
    fn <- function(n){rmvt(n, delta=u, sigma=sigmat, df=10, type="shifted")}
    return(fn)
}

get.d.centered.proposal <- function(u,sigmat){
    fn <- function(x){dmvt(x, delta=u, sigma=sigmat, df=10, log=F)}
    return(fn)
}

compute.weight.i <- function(i, xi, constrained.targets){
    #u <- optim(centers[i,], constrained.targets[[i]], control=list(fnscale=-1))$par
    u <- colMeans(xi)
    sigmat <- 3.5* diag(diag(cov(xi)))
    rq2 <- get.r.centered.proposal(u, sigmat)
    dq2 <- get.d.centered.proposal(u, sigmat)
    x2 <- rq2(25000)
    res.i <- bridge.sampling.very.fast(constrained.targets[[i]], dq2, xi, x2)
    return(list(res.i, x2))
}


# compute spectral clustering weights

w1.sc <- compute.weight.i(1, res.1.sc$X, sc.constrained.targets)[[1]]
w2.sc <- compute.weight.i(2, res.2.sc$X, sc.constrained.targets)[[1]]

ws.sc <- c(w1.sc, w2.sc) / (w1.sc + w2.sc)

save(ws.sc, sc, sc.indicator.fn, sc.constrained.targets, isshape1, isshape2, dfn, rprop.1, dprop.1,
     file="results/exploration-sc.rdata")


# compute voronoi clustering weights

w1.vc <- compute.weight.i(1, res.1.vc$X, vc.constrained.targets)[[1]]
w2.vc <- compute.weight.i(2, res.2.vc$X, vc.constrained.targets)[[1]]

ws.vc <- c(w1.vc, w2.vc) / (w1.vc + w2.vc)

save(vc.centers, ws.vc, vc, vc.indicator.fn, vc.constrained.targets, isshape1, isshape2, dfn, rprop.1, dprop.1,
     file="results/exploration-vc.rdata")


# # # # # # # # 
#
#res.sc.x <- rbind(res.1.sc$X, res.2.sc$X)
#plot(res.sc.x[,2] ~ res.sc.x[,1], xlim=c(-2,2), ylim=c(-2,2))
#
#x2.1 <- compute.weight.i(1, res.1.sc$X, sc.constrained.targets)[[2]]
#x2.2 <- compute.weight.i(2, res.2.sc$X, sc.constrained.targets)[[2]]
#x2.s <- rbind(x2.1, x2.2)
#
#plot(res.sc.x[,2] ~ res.sc.x[,1], xlim=c(-2,2), ylim=c(-2,2))
#points(x2.s[,1], x2.s[,2], col='red')
#
#res.vc.x <- rbind(res.1.vc$X, res.2.vc$X)
#plot(res.vc.x[,2] ~ res.vc.x[,1], xlim=c(-2,2), ylim=c(-2,2))
#
## # # # # # # #
#
#save(vc, indicator.fn, constrained.targets, isshape1, isshape2, dfn, rprop.1, dprop.1,
#     file="results/exploration-vc.rdata")
#
#clust.init <- get.clust.init(xs)
#clust.1.init <- clust.init[[1]]
#clust.2.init <- clust.init[[2]]
#
#
##plot(xs[,2] ~ xs[,1], col=vc$cluster)
#
#
#vc.constrained.targets <- get.constrained.targets(2, vc$indicator, dfn)
#res.vc.1 <- run.mv.mh(clust.1.init, 5000, rprop.1, dprop.1, vc.constrained.targets[[1]])
#res.vc.2 <- run.mv.mh(clust.2.init, 5000, rprop.1, dprop.1, vc.constrained.targets[[2]])
#
#
#
#
#res.sc.1 <- run.mv.mh(sc$centers[1,], 5000, rprop.1, dprop.1, constrained.targets[[1]])
#res.sc.2 <- run.mv.mh(sc$centers[2,], 5000, rprop.1, dprop.1, constrained.targets[[2]])
#
#res.sc.x <- rbind(res.sc.1$X, res.sc.2$X)
#plot(res.sc.x[,2] ~ res.sc.x[,1])
#
#plot(sc$Y[,2] ~ sc$Y[,1], col=sc$clust.clust)
#
#vc <- voronoi.clustering(xs, 2)
#
#plot(xs[,2] ~ xs[,1], col=vc$cluster)
#
#
#get.clust.init <- function(X){
#    clust.1.idx <- which(vc$cluster==1)
#    clust.2.idx <- which(vc$cluster==2)
#    return(list(X[sample(clust.1.idx,1),], X[sample(clust.2.idx,1),]))
#}
#
#clust.init <- get.clust.init(xs)
#clust.1.init <- clust.init[[1]]
#clust.2.init <- clust.init[[2]]
#
#vc.constrained.targets <- get.constrained.targets(2, vc$indicator, dfn)
#res.vc.1 <- run.mv.mh(clust.1.init, 5000, rprop.1, dprop.1, vc.constrained.targets[[1]])
#res.vc.2 <- run.mv.mh(clust.2.init, 5000, rprop.1, dprop.1, vc.constrained.targets[[2]])
#
#res.vc.x <- rbind(res.vc.1$X, res.vc.2$X)
#
#plot(res.vc.x[,2]~res.vc.x[,1])
#
#x1 <- do.call('rbind', lapply(seq(1000), function(i) rshape1()))
#x2 <- do.call('rbind', lapply(seq(1000), function(i) rshape2()))
#x <- rbind(x1,x2)
#plot(x[,2] ~ x[,1], xlim=c(-2,2), ylim=c(-2,2), col=res$cluster)
#
