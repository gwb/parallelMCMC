

importance.sampling <- function(dtarget, r.importance.distr, d.importance.distr, niter){
  X <- r.importance.distr(niter)
  return(1/niter * sum(apply(X, 1, function(x) dtarget(x)/ d.importance.distr(x))))
}


geometric.bridge <- function(dq1, dq2){
  fn <- function(x){
    if(dq1(x) == 0 || dq2(x) == 0){
      return(0)
    }else{
      return( 1 / sqrt(dq1(x) * dq2(x)))
    }
  }
  return(fn)
}

bridge.sampling <- function(dq1, dq2, x1, x2, bridge.fn){
  numerator <- 1 / length(x2) * sum(apply(x2, 1, function(x) dq1(x)*bridge.fn(x)))
  denominator <- 1 / length(x1) * sum(apply(x1,1, function(x) dq2(x)*bridge.fn(x)))
  return(numerator/denominator)
}


####

