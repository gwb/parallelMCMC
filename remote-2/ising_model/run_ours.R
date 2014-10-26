source('functions.R')

### Loading exploration results
load("results/exploration.rdata")

args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[1])

### Begin Parameters
N <- 3
Beta <- 1 / 1.3
J <- 1
analytical.distr <- c(0.4892646784, 0.0093584720, 0.0009026467, 0.0003708764, 0.0001033265,
                      0.0001033265, 0.0003708764,  0.0009026467, 0.0093584720, 0.4892646784)

### End Parameter


### Run Parallel Chains
print("Running parallel chains")

res.mh.1 <- run.mh.modif(M.1, prop.spin, get.A.i[[1]], Beta=Beta, niter=10000)
res.magnet.1 <- get.magnetization.distr(res.mh.1[[1]])

res.mh.2 <- run.mh.modif(M.2, prop.spin, get.A.i[[2]], Beta=Beta, niter=10000)
res.magnet.2 <- get.magnetization.distr(res.mh.2[[1]])

res.magnet <- (ws[1] * res.magnet.1[,2] + ws[2] * res.magnet.2[,2]) / (ws[1] * sum(res.magnet.1) + ws[2] * sum(res.magnet.2))

hellinger.distance <- (1/ sqrt(2)) * sqrt(sum((sqrt(res.magnet) - sqrt(analytical.distr))^2))

save(hellinger.distance, file=paste("results/ours/ours",idx, ".rdata", sep=""))

