

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

