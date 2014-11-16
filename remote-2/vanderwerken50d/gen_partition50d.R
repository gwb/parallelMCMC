source("model50.R")
source("../../mh.R", chdir=T)
source("../../graphics.R", chdir=T)
source("../../bridge-sampling.R", chdir=T)
source("../../parallel-tempering.R", chdir=T)

set.seed(123)

sd.sig.prop <- sqrt(0.06) # standard deviation of the proposal
sd.sig.prop.2 <- sqrt(0.12)
sd.sig.prop.3 <- sqrt(0.2)
sd.sig.prop.4 <- sqrt(0.4)
sd.sig.prop.5 <- sqrt(0.7)

dprop <- function(nstate, ostate){
    return(prod(dnorm(nstate, ostate, sd.sig.prop)))
}

dprop.2 <- function(nstate, ostate){
    return(prod(dnorm(nstate, ostate, sd.sig.prop.2)))
}

dprop.3 <- function(nstate, ostate){
    return(prod(dnorm(nstate, ostate, sd.sig.prop.3)))
}

dprop.4 <- function(nstate, ostate){
    return(prod(dnorm(nstate, ostate, sd.sig.prop.4)))
}

dprop.5 <- function(nstate, ostate){
    return(prod(dnorm(nstate, ostate, sd.sig.prop.5)))
}

rprop <- function(ostate){
    return(rnorm(length(ostate), mean=ostate,sd=sd.sig.prop))
}

rprop.2 <- function(ostate){
    return(rnorm(length(ostate), mean=ostate,sd=sd.sig.prop.2))
}

rprop.3 <- function(ostate){
    return(rnorm(length(ostate), mean=ostate,sd=sd.sig.prop.3))
}

rprop.4 <- function(ostate){
    return(rnorm(length(ostate), mean=ostate,sd=sd.sig.prop.4))
}

rprop.5 <- function(ostate){
    return(rnorm(length(ostate), mean=ostate,sd=sd.sig.prop.5))
}

dtarget.2 <- function(x) dtarget(x)^(0.75)
dtarget.3 <- function(x) dtarget(x)^(0.5)
dtarget.4 <- function(x) dtarget(x)^(0.25)
dtarget.5 <- function(x) dtarget(x)^(0.12)


single.dproposals.ls <- list(dprop, dprop.2, dprop.3, dprop.4, dprop.5)
single.rproposals.ls <- list(rprop, rprop.2, rprop.3, rprop.4, rprop.5)
dtarget.ls <- list(dtarget, dtarget.2, dtarget.3, dtarget.4, dtarget.5)

X0 <- do.call('rbind', lapply(seq(length(dtarget.ls)), function(i) u.ls[[1]]))

res <- parallel.tempering.fast(dtarget.ls, single.rproposals.ls, single.dproposals.ls, X0, niter=25000)

states <- res[[1]][[1]]
#algo.dt <- plot.vhd.no.clusters(states,nproj=10)
#ggplot(data=algo.dt[[1]], aes(x=x,y=y)) + geom_point() + facet_grid(.~proj)

Y <- states[sample(seq(nrow(states)), 10000),]

vc <- kmeans(Y, 4, iter.max=100, nstart=20)

indicator.fn <- function(x){
    return(as.vector(which.min(apply(vc$centers, 1, function(y) sum( (x-y)^2 )))))
}


# our method

constrained.targets <- get.constrained.targets(4, indicator.fn, dtarget)


get.r.centered.proposal <- function(u, sigmat){
    #fn <- function(n){rmvt(n, delta=u, sigma=sigmat, df=10, type="shifted")}
    fn <- function(n) {rmvnorm(n, mean=u, sigma=1.5*sigmat)}
    return(fn)
}

get.d.centered.proposal <- function(u,sigmat){
    #fn <- function(x){dmvt(x, delta=u, sigma=sigmat, df=10, log=F)}
    fn <- function(x){dmvnorm(x, mean=u, sigma=1.5*sigmat, log=F)}
    return(fn)
}

compute.weight.i <- function(i, xi){
    #u <- optim(vc$centers[i,], constrained.targets[[i]], control=list(fnscale=-1))$par
    u <- vc$centers[i,]
    sigmat <- cov(xi)
    #sigmat <- 0.5*diag(length(vc$centers[i,]))
    rq2 <- get.r.centered.proposal(u, sigmat)
    dq2 <- get.d.centered.proposal(u, sigmat)
    x2 <- rq2(200000)
    res.i <- bridge.sampling.very.fast(constrained.targets[[i]], dq2, xi, x2)
    return(list(res.i))
}

c.run.simul <- function(){
    res.1 <- run.mv.mh.fast(vc$centers[1,], 30000, rprop, dprop, constrained.targets[[1]])
    res.2 <- run.mv.mh.fast(vc$centers[2,], 30000, rprop, dprop, constrained.targets[[2]])
    res.3 <- run.mv.mh.fast(vc$centers[3,], 30000, rprop, dprop, constrained.targets[[3]])
    res.4 <- run.mv.mh.fast(vc$centers[4,], 30000, rprop, dprop, constrained.targets[[4]])
    res.ls <- list(res.1$X, res.2$X, res.3$X, res.4$X)
    wis <- unlist(lapply(seq(4), function(i) compute.weight.i(i, res.ls[[i]])))
    wis <- wis / sum(wis)
    return(wis)
}


ws <- c.run.simul()

save(C, ws, dtarget, dprop, rprop, constrained.targets, vc, u.ls, sd.sig.prop, dmvn, indicator.fn, Y,file="results/exploration.rdata")

