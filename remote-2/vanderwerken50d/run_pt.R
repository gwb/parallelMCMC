require(methods)
source("../../mh.R", chdir=T)
source("../../parallel-tempering.R", chdir=T)

# # # # # # # # # # # # # # # # 

C <- 4
d <- 50

load("results/uls.rdata")

dmvn <- function(x,i){
    return(prod(dnorm(x, u.ls[[i]], sqrt(0.5))))
}

dtarget <- function(x){
    w.ls <- c(0.1, 0.2, 0.3, 0.4)
    res <- 0
    for(i in seq(C)){
        res <- res + w.ls[i] * dmvn(x, i)
    }
    return(res)
}

# # # # # # # # # # # # # # # # 

args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[1])



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

res <- parallel.tempering.fast(dtarget.ls, single.rproposals.ls, single.dproposals.ls, X0, niter=65000)

states <- res[[1]][[1]]


u.hat.pt <- colMeans(states)


save(u.hat.pt, file=paste("results/pt/pt", idx, ".rdata", sep=""))
