source("../models/vanderWerken-1.R")
source("parallel_adaptive-hd.R")
source("../bridge-sampling.R")
source("../graphics.R")
source("../clustering.R")


# # # # # # # # # # # #
# Running the sampler
# # # # # # # # # # # #

# Number of partition-sets we want (should correspond to number of parallel cores we have)
K <- 4

# Initial centers for exploration phase

init.points <- matrix(c(7, -3,
                        2,7), byrow=T, nrow=2)



mh.explorer <- get.mv.mh(2500, draw.normal.proposal, eval.normal.proposal, burnin=1000)

foo <- init.exploration.test(dtarget, mh.explorer, init.points)

bar <- init.exploration(dtarget, mh.explorer, init.points)

system.time(init.exploration.test(dtarget, mh.explorer, init.points))

system.time(init.exploration(dtarget, mh.explorer, init.points))

system.time(init.exploration.old(dtarget, mh.explorer, matrix(init.points[1,],nrow=1)))


init.points <- matrix(c(-4, -4,
                        -4,0,
                        -4,4,
                        0,-4,
                        0,0,
                        0,4,
                        4,-4,
                        4,0,
                        4,4), byrow=T, nrow=9)

mh.explorer <- get.mv.mh(1000, draw.normal.proposal, eval.normal.proposal, burnin=100)

system.time(init.exploration.test(dtarget, mh.explorer, init.points))

system.time(init.exploration(dtarget, mh.explorer, init.points))


