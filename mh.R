run.mh <- function(x0, niter, rproposal, dproposal, target){
  X <- c()
  L <- c()
  R <- c()
  X.tm1 <- x0
  for(i in seq(niter)){
    L.t <- rproposal(X.tm1)
    r <- (target(L.t) * dproposal(X.tm1, L.t)) / (target(X.tm1) * dproposal(L.t, X.tm1))

    X <- c(X, X.tm1)
    
    if(runif(1) < r){
      #X <- c(X,L.t)
      X.tm1 <- L.t
    } #else {
      #X <- c(X, X.tm1)
    #}
    
    R <- c(R,r)
    L <- c(L, L.t)
  }

  acceptance.rate <- length(unique(X))/ length(X)
  
  return(list(cbind(X, L, R), acceptance.rate))
}
