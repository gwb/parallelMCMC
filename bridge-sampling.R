

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



geometric.bridge.fast <- function(dq1, dq2){
  fn <- function(x.ls){
    dq2x.ls <- dq2(x.ls)
    dq1x.ls <- apply(x.ls, 1, dq1)
    return(ifelse(dq2x.ls * dq1x.ls==0, 0, 1/sqrt(dq1x.ls*dq2x.ls)))
  }
  return(fn)
}

bridge.sampling.fast <- function(dq1, dq2, x1, x2, bridge.fn){
  dq1x2.ls <- apply(x2, 1,dq1)#dq1(x2)
  dq2x1.ls <- dq2(x1)
  bridge.x1.ls <- bridge.fn(x1)
  bridge.x2.ls <- bridge.fn(x2)
  numerator <- 1/length(x2) * sum(dq1x2.ls * bridge.x2.ls)
  denominator <- 1/length(x1) * sum(dq2x1.ls * bridge.x1.ls)
  return(numerator/denominator)
}


bridge.sampling.very.fast <- function(dq1, dq2, x1, x2){
  dq1x2.ls <- apply(x2, 1, dq1)#dq1(x2)
  dq1x1.ls <- apply(x1, 1, dq1)#dq1(x1)
  dq2x2.ls <- dq2(x2)
  dq2x1.ls <- dq2(x1)
  # geometric bridges
  bridge.x1.ls <- ifelse(dq2x1.ls * dq1x1.ls==0, 0, 1/sqrt(dq2x1.ls * dq1x1.ls))
  bridge.x2.ls <- ifelse(dq2x2.ls * dq1x2.ls==0, 0, 1/sqrt(dq2x2.ls * dq1x2.ls))
  numerator <- 1/length(x2) * sum(dq1x2.ls * bridge.x2.ls)
  denominator <- 1/length(x1) * sum(dq2x1.ls * bridge.x1.ls)
  return(numerator/denominator)
}

####

