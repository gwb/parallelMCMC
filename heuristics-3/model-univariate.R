

# note: the model we consider is for:
#        f(x) = 0.5 * N( x | -mu, 1) + 0.5 * N( x | mu, 1)

get.dtarget <- function(mu){
    dtarget <- function(x){
        return( 0.5 * dnorm(x, -mu, 1) + 0.5 * dnorm(x, mu, 1) )
    }
    return(dtarget)
}

get.dproposal <- function(tau){
    dproposal <- function(x, x.star){
        return( 1/(2*tau) * (abs(x-x.star) < tau))
    }
    return(dproposal)
}

get.dtransition <- function(tau, mu, dprop, dtarget){
    dtransition <- function(x, x.star){
        return( dprop(x, x.star) * min(1, dtarget(x.star) / dtarget(x)) )
    }
}

rtarget <- function(n, mu){
    res <- vector('numeric', n)
    for(i in seq(n)){
        if(runif(1) < 0.5){
            res[i] <- rnorm(1, -mu, 1)            
        } else {
            res[i] <- rnorm(1, mu, 1)
        }
    }
    return(res)
}

get.rproposal <- function(tau){
    rproposal <- function(x){
        return(x + runif(1, -tau, tau))
    }
    return(rproposal)
}

get.weights <- function(mu, mid){
    p.mid.val <- 0.5 * pnorm(mid, -mu, 1) + 0.5 * pnorm(mid, mu, 1)
    return(c(p.mid.val, 1-p.mid.val))
}
