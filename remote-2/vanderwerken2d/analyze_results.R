


load('results/exploration.rdata')

true.mean <- 0.02 * u.ls[[1]] + 0.2 * u.ls[[2]] + 0.2 * u.ls[[3]] + 0.58 * u.ls[[4]]


# Our method

sim.filelist <- list.files('results/vw/')

mu.hat.ls <- NULL
for(filename in sim.filelist){
    load(file.path('results/vw', filename))
    mu.hat.ls <- rbind(mu.hat.ls, mu.hat)
}

dimnames(mu.hat.ls) <- NULL

dist.sq <- apply(mu.hat.ls, 1, function(x) sum((x-true.mean)^2))

mean.dist.sq <- mean(dist.sq)
var.dist.sq <- var(dist.sq)


save(mu.hat.ls, true.mean, dist.sq, mean.dist.sq, var.dist.sq, file='results/vw_analysis.data')


# parallel tempering method

pt_directory <- 'results/vw_pt'
sim.filelist <- list.files(pt_directory)
mu.hat.ls <- NULL

for(filename in sim.filelist){
    load(file.path(pt_directory, filename))
    mu.hat.ls <- rbind(mu.hat.ls, mu.hat)
}

dimnames(mu.hat.ls) <- NULL

dist.sq <- apply(mu.hat.ls, 1, function(x) sum((x-true.mean)^2))

mean.dist.sq <- mean(dist.sq)
var.dist.sq <- var(dist.sq)

save(mu.hat.ls, true.mean, dist.sq, mean.dist.sq, var.dist.sq, file='results/pt_analysis.data')


