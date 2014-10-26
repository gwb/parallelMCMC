source('functions.R')


args <- commandArgs(trailingOnly=TRUE)
idx <- as.integer(args[1])

### Begin Parameters
N <- 3
Beta <- 1 / 1.3
J <- 1
analytical.distr <- c(0.4892646784, 0.0093584720, 0.0009026467, 0.0003708764, 0.0001033265,
                      0.0001033265, 0.0003708764,  0.0009026467, 0.0093584720, 0.4892646784)
### End Parameter


### Begin Parallel Tempering
print("Running parallel tempering")
Beta.ls <- 1 / c(Beta, 2, 2.4)
M.0 <- matrix(rep(1,N^2), nrow=N)#matrix(sample(c(-1,1), size=N^2, replace=T), nrow=N)
M.0.ls <- lapply(seq(length(Beta.ls)), function(i) M.0)
res.pt <- parallel.tempering(M.0.ls, prop.spin, .get.A.i, get.swap.A.i, Beta.ls, niter=40000, N.spin = N)
### End Parallel Tempering


M.ls <- lapply(seq(nrow(res.pt[[1]][[1]])), function(i) matrix(res.pt[[1]][[1]][i,], nrow=N))
magnet <- get.magnetization.distr(M.ls)
res.magnet <- magnet[,2] / sum(magnet[,2])

hellinger.distance <- (1/ sqrt(2)) * sqrt(sum((sqrt(res.magnet) - sqrt(analytical.distr))^2))


save(hellinger.distance, file=paste("results/pt/pt",idx, ".rdata", sep=""))
