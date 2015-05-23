source('functions.R')
source('../../bridge-sampling.R', chdir=T)

### Begin Parameters
N <- 3
kbT <- 2.2
Beta <- 1 / kbT
J <- 1
### End Parameter


### Begin Parallel Tempering
print("Running parallel tempering")
Beta.ls <- 1 / c(kbT, 2.25, 2.4)
M.0 <- matrix(rep(1,N^2), nrow=N)#matrix(sample(c(-1,1), size=N^2, replace=T), nrow=N)
M.0.ls <- lapply(seq(length(Beta.ls)), function(i) M.0)
res.pt <- parallel.tempering(M.0.ls, prop.spin, .get.A.i, get.swap.A.i, Beta.ls, niter=30000, N.spin = N)
### End Parallel Tempering

### Begin Debug Parallel Tempering
#M.ls <- lapply(seq(nrow(res.pt[[1]][[1]])), function(i) matrix(res.pt[[1]][[1]][i,], nrow=N))
#magnet <- get.magnetization.distr(M.ls)
#barplot(magnet[,2], names=magnet[,1])
### End Debug Parallel Tempering

### Begin Spectral Clustering
print("Performing spectral clustering")
X <- tail(res.pt[[1]][[1]], 10000)
sc <- spectral.clustering.vec(X, nclust=2, Q, nsub=150)
indicator.fn <- function(x) sc$indicator(as.vector(x))
### End Spectral Clustering

### Begin Debug
#magnet.1 <- unlist(sapply(which(sc$cluster==1), function(i) sum(sc$Y[i,]) / N^2))
#magnet.2 <- unlist(sapply(which(sc$cluster==2), function(i) sum(sc$Y[i,]) / N^2))
#magnet.all <- c(magnet.1, magnet.2)
#plot(magnet.all~seq(length(magnet.all)), col=c(rep(1, length(magnet.1)), rep(2, length(magnet.2))))
#### End Debug

print("Running parallel chains")

get.A.i <- get.constrained.Ai(2, indicator.fn, .get.A.i)

M.1 <- matrix(rep(1, N^2), nrow=N)
M.2 <- -M.1
if(get.A.i[[1]](M.1, c(1,1), 0.5) == 0){
    M.2 <- M.1
    M.1 <- -M.1
}

res.mh.1 <- run.mh.modif(M.1, prop.spin, get.A.i[[1]], Beta=Beta, niter=30000)
res.magnet.1 <- get.magnetization.distr(res.mh.1[[1]])
#barplot(res.magnet.1[,2], names=res.magnet.1[,1])

res.mh.2 <- run.mh.modif(M.2, prop.spin, get.A.i[[2]], Beta=Beta, niter=30000)
res.magnet.2 <- get.magnetization.distr(res.mh.2[[1]])
#barplot(res.magnet.2[,2], names=res.magnet.2[,1])

### End Parallel

### Begin Bridge Sampling
print("Performing bridge sampling")
props.1 <- get.proposals(res.mh.1, cutoff=0.995)
r.prop.1 <- props.1[[1]]
d.prop.1 <- props.1[[2]]

props.2 <- get.proposals(res.mh.2)
r.prop.2 <- props.2[[1]]
d.prop.2 <- props.2[[2]]


n.unnorm.dens.vec <- function(x) unnorm.dens.vec(x, Beta)

x1.1 <- do.call('rbind', lapply(res.mh.1[[1]], function(x) as.vector(x)))
x2.1 <- do.call('rbind', lapply(seq(3000), function(i) r.prop.1()))
alph.1 <- geometric.bridge(n.unnorm.dens.vec, d.prop.1)
foo.1 <- bridge.sampling(n.unnorm.dens.vec, d.prop.1, x1.1, x2.1, alph.1)

x1.2 <- do.call('rbind', lapply(res.mh.2[[1]], function(x) as.vector(x)))
x2.2 <- do.call('rbind', lapply(seq(3000), function(i) r.prop.2()))
alph.2 <- geometric.bridge(n.unnorm.dens.vec, d.prop.2)
foo.2 <- bridge.sampling(n.unnorm.dens.vec, d.prop.2, x1.2, x2.2, alph.2)


ws <- c(foo.1, foo.2) / (foo.1 + foo.2)
### End Bridge Sampling

save(ws, sc, indicator.fn, get.A.i, res.mh.1, res.magnet.1, res.mh.2, res.magnet.2, M.1, M.2, Beta, file="results/exploration.rdata")

