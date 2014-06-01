
require(ggplot2)

# # # # # # #
# FUNCTIONS #
# # # # # # #

get.numerical.intervals <- function(intervals){
  res <- NULL
  for(interval in intervals){
    res <- rbind(res, as.numeric(strsplit(substr(interval, 2, nchar(interval) - 1), ',')[[1]]))
  }
  return(res)
}

add.rects <- function(clust.ints, my.colors){
  lim <- par('usr')
  for(i in seq(nrow(clust.ints))){
    rect(clust.ints[i,1], lim[3]-1, clust.ints[i,2], lim[4]+1, col=my.colors[clust.ints[i,3]], border=NA)
  }
}

# color scheme with transparency
my.cols <- c(rgb(1,0,0,alpha=0.2), rgb(0,1,0, alpha=0.2), rgb(0,0,1, alpha=0.2))


# plot partition (for 2D data)
get.partition.colors <- function(xs, ys, cluster.fn){
  res <- expand.grid(xs, ys)
  res <- cbind(res, apply(res, 1, cluster.fn))
  names(res) <- NULL
  return(res)
}


# plots the contours of a bivariate density 
plot.bivariate.density <- function(fn, x1range=c(-10,10), x2range=c(-10,10), bins.per.dim=100){
  x1 <- seq(x1range[1], x1range[2], length.out = bins.per.dim)
  x2 <- seq(x2range[1], x2range[2], length.out = bins.per.dim)
  
  z2 <- matrix(0, nrow=bins.per.dim, ncol=bins.per.dim)
  
  for(i in seq(1,bins.per.dim)){
    for(j in seq(1,bins.per.dim)){
      z2[i,j] <- fn(c(x1[i], x2[j]))
    }
  }
  
  contour(list(x=x1,y=x2,z=z2), col='red',add=F)

}

plot.cluster.and.target <- function(cluster.fn, rtarget, x1range=c(-10,10), x2range=c(-10,10), target.draws=10000){
  pcolors <- get.partition.colors(seq(x1range[1],x1range[2]), seq(x2range[1], x2range[2]), cluster.fn)
  dt <- data.frame(pcolors)
  names(dt) <- c('x','y','z')

  dt2 <- data.frame(rtarget(target.draws))
  names(dt2) <- c('x','y')

  g <- ggplot(data=dt2, aes(x,y)) + stat_density2d() + geom_tile(data=dt, aes(x,y,fill=factor(z), alpha=0.2))
  return(g)
}

# # # # # # #
# EXAMPLES  #
# # # # # # #

# The examples below show how to plot some quantities of interest

# -> the target
# curve(dtarget, from=min.val, to = max.val, n=5000)

# -> kernel density estimation, from the draws
# plot(density(all.res[,1], n=2000, bw="bcv"))

# ->superimposing the target and the histogram
# curve(dtarget, from=min.val, to = max.val, n=5000, add=FALSE, col="blue")
# hist(all.res[,1], breaks=1000, col="red", freq=FALSE, add=TRUE)
# curve(dtarget, from=min.val, to = max.val, n=5000, add=TRUE, col="blue", lwd=2)

# -> progressively adding histogram (that's pretty cool)
#     (proof of concept for live monitoring of mcmc)
#
# rand.res = sample(all.res[,1], length(all.draws))
# for(i in seq(1, 50000, by=100)){
#   curve(dtarget, from=0, to = 10, n=5000, add=FALSE, col="blue")
#   hist(c(rand.res[seq(1, i)], rep(10,50000-i)), breaks=1000, col="red", freq=FALSE, add=TRUE)
#  S ys.sleep(0.5)
# }

#
# cfn <- proj.clustering(X, 1, 2)
# bli <- get.partition.colors(seq(-1, 5), seq(-1,5), cfn)
# dt <- data.frame(bli)
# names(dt) <- c('x','y','z')
# qplot(x,y,fill=factor(z), data=dt, geom="tile")


#
#
#
#
#
#
