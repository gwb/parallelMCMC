require(methods)
require(mvtnorm)
source("../../parallel-tempering.R", chdir=T)
source("../../mh.R", chdir=T)
source("../../clustering.R", chdir=T)
source("../../bridge-sampling.R", chdir=T)
source("../../graphics.R", chdir=T)
require(mvtnorm)

args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[1])

load('results/exploration.rdata')

# Define target and all

C <- 4
d <- 10

#sig.ls <- vector('list', C)
#u.ls <- vector('list', C)

w.ls <- c(0.1, 0.2, 0.3, 0.4)

#sig.prop <- 0.2*diag(d)

#for(i in seq(C)){
#    sig.ls[[i]] <- 0.5*diag(d)
#    u.ls[[i]] <- runif(d, min=-10, max=10)
#}


# targets

dtarget <- function(x){
      res <- 0
        for(i in seq(C)){
                res <- res + w.ls[i] * dmvnorm(x, u.ls[[i]], sig.ls[[i]])
            }
        return(res)
  }

rtarget <- function(n){
      res <- NULL
        for(j in seq(n)){
                i <- sample(seq(C), 1, prob=w.ls)
                res <- rbind(res, rmvnorm(1, u.ls[[i]], sig.ls[[i]]))
            }
        return(res)
  }


dfn.1 <- function(x) {dtarget(x)}
dfn.2 <- function(x) {dtarget(x)^(0.75)}
dfn.3 <- function(x) {dtarget(x)^(0.5)}
dfn.4 <- function(x) {dtarget(x)^(0.1)}
dfn.5 <- function(x) {dtarget(x)^(0.05)}
dfn.6 <- function(x) {dtarget(x)^(0.01)}


dprop.1 <- function(nstate, ostate){
    return(dmvnorm(nstate, mean=ostate, sigma=sig.prop))
}
dprop.2 <- function(nstate, ostate){
        return(dmvnorm(nstate, ostate, 2 * sig.prop))
    }
dprop.3 <- function(nstate, ostate){
        return(dmvnorm(nstate, ostate, 5 * sig.prop))
    }
dprop.4 <- function(nstate, ostate){
        return(dmvnorm(nstate, ostate, 7 * sig.prop))
    }
dprop.5 <- function(nstate, ostate){
        return(dmvnorm(nstate, ostate, 10 * sig.prop))
    }
dprop.6 <- function(nstate, ostate){
        return(dmvnorm(nstate, ostate, 14 * sig.prop))
    }



rprop.1 <- function(ostate){
    return(rmvnorm(1, mean=ostate, sigma=sig.prop))
}
rprop.2 <- function(ostate){
        return(rmvnorm(1, ostate, 2 * sig.prop))
    }
rprop.3 <- function(ostate){
        return(rmvnorm(1, ostate, 5 * sig.prop))
    }
rprop.4 <- function(ostate){
        return(rmvnorm(1, ostate, 7 * sig.prop))
    }
rprop.5 <- function(ostate){
        return(rmvnorm(1, ostate, 10 * sig.prop))
    }
rprop.6 <- function(ostate){
        return(rmvnorm(1, ostate, 14 * sig.prop))
    }


dprop.1.vec <- function(nstate, ostate){
    if(!is.null(dim(ostate))){
        if(dim(ostate)[1] == 1){
            ostate <- as.vector(ostate)
        }
        return(dprop.1(ostate, nstate))
    } else if(!is.null(dim(nstate))){
        if(dim(nstate)[1] == 1){
            nstate <- as.vector(nstate)
        }
        return(dprop.1(nstate, ostate))
    } else {
        return(dprop.1(nstate, ostate))
    }
}



# Exploring with parallel tempering
targets.ls <- list(dfn.1, dfn.2, dfn.3, dfn.4, dfn.5, dfn.6)
single.rproposals.ls <- list(rprop.1, rprop.2, rprop.3, rprop.4, rprop.5, rprop.6)
single.dproposals.ls <- list(dprop.1, dprop.2, dprop.3, dprop.4, dprop.5, dprop.6)
init.states.mat <- matrix(rep(u.ls[[1]], 6), nrow=6, byrow=T)

#set.seed(123)

#res <- parallel.tempering(targets.ls, single.rproposals.ls, single.dproposals.ls, init.states.mat, 20000)

res <- parallel.tempering.notrace(targets.ls, single.rproposals.ls, single.dproposals.ls, init.states.mat, 35000)
states <- res[[1]]

mu.hat <- colMeans(states)

save(mu.hat, file=paste("results/pt/pt", idx, ".rdata", sep=""))

