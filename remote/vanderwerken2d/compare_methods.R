
load("results/good.rdata")
good.res <- res

load("results/bad.rdata")
bad.res <- res

load("results/exploration.rdata")

true.mean <- ws[1] * u.ls[[1]] + ws[2] * u.ls[[2]] + ws[3] * u.ls[[3]] + ws[4] * u.ls[[4]]

colMeans(good.res)
mean(apply(good.res, 1, function(x) sum( (x-true.mean)^2 )))

colMeans(bad.res)
mean(apply(bad.res, 1, function(x) sum( (x-true.mean)^2 )))
