
# all simultions where run for a total of 5000 iterations 

require(ggplot2)

load("results/good.rdata")
good.res <- res

load("results/bad.rdata")
bad.res <- res

load("results/hmc.rdata")
hmc.res <- res

load("results/exploration.rdata")

true.mean <- ws[1] * u.ls[[1]] + ws[2] * u.ls[[2]] + ws[3] * u.ls[[3]] + ws[4] * u.ls[[4]]

colMeans(good.res)
good.dist <- apply(good.res, 1, function(x) sum( (x-true.mean)^2 ))

colMeans(bad.res)
bad.dist <- apply(bad.res, 1, function(x) sum( (x-true.mean)^2 ))

colMeans(hmc.res)
hmc.dist <- apply(hmc.res, 1, function(x) sum( (x-true.mean)^2 ))


sds <- c(sd(good.dist), sd(bad.dist), sd(hmc.dist))
ms <- c(mean(good.dist), mean(bad.dist), mean(hmc.dist))
emax <- ms + sds
emin <- ms - sds
ys <- factor(c("ours", "parallel tempering", "stan(hmc)"))

dt <- data.frame(ms=ms,
                 emax=emax,
                 emin=emin,
                 ys=ys)

ggplot(data=dt, aes(ms, ys)) + geom_point() + geom_errorbarh(aes(xmax=emax, xmin=emin, height=0.1))



