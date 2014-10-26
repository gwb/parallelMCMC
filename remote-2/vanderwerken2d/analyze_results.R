


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

vw.dist.sq <- apply(mu.hat.ls, 1, function(x) sum((x-true.mean)^2))

vw.mean.dist.sq <- mean(vw.dist.sq)
vw.var.dist.sq <- var(vw.dist.sq)


save(mu.hat.ls, true.mean, vw.dist.sq, vw.mean.dist.sq, vw.var.dist.sq, file='results/vw_analysis.data')

print("Our method:")
print(paste("----- mean: ", vw.mean.dist.sq))
print(paste("----- sd:   ", sqrt(vw.var.dist.sq)))

# parallel tempering method

pt_directory <- 'results/vw_pt/'
sim.filelist <- list.files(pt_directory)
mu.hat.ls <- NULL

for(filename in sim.filelist){
    load(file.path(pt_directory, filename))
    mu.hat.ls <- rbind(mu.hat.ls, mu.hat)
}

dimnames(mu.hat.ls) <- NULL

pt.dist.sq <- apply(mu.hat.ls, 1, function(x) sum((x-true.mean)^2))

pt.mean.dist.sq <- mean(pt.dist.sq)
pt.var.dist.sq <- var(pt.dist.sq)

save(mu.hat.ls, true.mean, pt.dist.sq, pt.mean.dist.sq, pt.var.dist.sq, file='results/pt_analysis.data')

print("Parallel Tempering:")
print(paste("----- mean: ", pt.mean.dist.sq))
print(paste("----- sd:   ", sqrt(pt.var.dist.sq)))

# Naive method

naive <- directory <- 'results/naive/'
sim.filelist <- list.files(naive <- directory)
mu.hat.ls <- NULL

for(filename in sim.filelist){
        load(file.path(naive <- directory, filename))
            mu.hat.ls <- rbind(mu.hat.ls, mu.hat)
    }

dimnames(mu.hat.ls) <- NULL

naive.dist.sq <- apply(mu.hat.ls, 1, function(x) sum((x-true.mean)^2))
naive.mean.dist.sq <- mean(naive.dist.sq)
naive.var.dist.sq <- var(naive.dist.sq)

save(mu.hat.ls, true.mean, naive.dist.sq, naive.mean.dist.sq, naive.var.dist.sq, file='results/naive_analysis.data')


print("Naive parallel:")
print(paste("----- mean: ", naive.mean.dist.sq))
print(paste("----- sd:   ", sqrt(naive.var.dist.sq)))
