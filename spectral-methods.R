

# # # # # # #
# FUNCTIONS #
# # # # # # #


#
# Allows us to test different types of partitioning algorithms. If A is asym, consider:
# 1) t(A) %*% A
# 2) t(A) + A
# Both of which are symmetric, and support the following 
#
.sym.spectral.clustering <- function(G, K){
  D.minus.sqrt <- diag(1/sqrt(rowSums(G)))
  L <- diag(nrow(G)) - D.minus.sqrt %*% G %*% D.minus.sqrt
  X <- eigen(L)$vectors[,tail(seq(nrow(L)), K)]
  T <- t(apply(X, 1, function(x) x/sqrt(sum(x^2))))
  return(kmeans(T,K,iter.max=50,nstart=20)$cluster)
}

# Applies spectral clustering to the symmetrized matrix t(G) %*% G
sym.spectral.clustering <- function(G,K){
  return(.sym.spectral.clustering(t(G) %*% G, K))
}

# The functions below are implementations of methods from Meila & Pentney
# asym.spectral.clustering - corresponds to the *BestWCut* algorithm
# wcut - corresponds to the quantity *WCut* defined in the paper
#

asym.spectral.clustering <- function(G, K){
  # following Meila & Pentney: Clustering by weighted cuts in directed graphs
  #
  # Note: the "eigen" function un R automatically returns orthonormal eigenvectors
  # when the input matrix is symmetric, so we don't have any additional step
  
  HB <- 0.5 * (2 * diag(nrow(G)) - G - t(G))
  X <- eigen(HB)$vectors[,tail(seq(nrow(HB)), K)]
  T <- t(apply(X, 1, function(x) x/sqrt(sum(x^2))))
  return(kmeans(T,K,iter.max=30,nstart=10)$cluster)
}


# The purpose of this function is to reassign the cluster-blocks that contain
# no information to the closest real clusters. During the clustering phase, there
# are some intervals in the discretization on which we have no information (they are
# not visited by any of the parallel chains), so the spectral clustering assigns them
# pretty much randomly, or not very meaningfully. This postprocessing function identifies
# those intervals, and reassigns them to the closest real cluster.
# Args - clusters: a vector of number in {1..K} indication for each position, which cluster
#                  they belong to
#      - K.mat   : Some informtations obtained from the second output of get.K.emp.by.index.new
#      - K       : the number of clusters 
post.process.clusters <- function(clusters, K.mat, K){
  tK.mat <- t(K.mat)
  noinfo <- which(tK.mat[,1] == 0 & tK.mat[,2] == 0)
  info <- which(tK.mat[,1] != 0 | tK.mat[,2] != 0)
  
  clusters.centers <- NULL
  info.clusters <- clusters[info]

  cdist <- function(x, c.center){
    return(abs(x - c.center))
  }
  
  for(i in seq(1,K)){
    clusters.i <- which(clusters == i)
    clusters.i.info <- clusters.i[clusters.i %in% info]
    clusters.centers = c(clusters.centers, mean(clusters.i.info))
  }

  for(j in noinfo){
    clusters[j] <- which.min(cdist(j, clusters.centers))
  }
  return(clusters)
}


# computes the value of the weighted cut for a certain cluster. Useful
# to check wether the *mincut* found by the algorithm is good or not. See
# Meila & Pentney for details about T and Tp
wcut <- function(A, T, Tp, clusters){
  # A is a (possibly) asymmetric adjacency matrix
  # T is a vector (of weight)
  # Tp is also a vector of weights
  # clusters is a vector where each entry indicates the cluser
  #          to which a node belongs
  # The function computes the WNCut as defined in Meila & Pentney

  res <- 0

  for(k in unique(clusters)){
    for(kp in unique(clusters)){
      if(k != kp){
        cut = 0
        for(i in which(clusters == k)){
          for(j in which(clusters == kp)){
            cut <- cut + Tp[i] * A[i,j]
          }
        }
        res <- res + cut / sum(T[which(clusters==k)])
      }
    }
  }
  return(res)
}


# # # # # # #
# EXAMPLES  #
# # # # # # #

# -> simple symmetric transition matrix, which should be clustered
#    as (1,2), and (3,4) if we allow 2 clusters
#
# test.m <- matrix(c(0.5 , 0.3 , 0.1 , 0.1,
#                    0.3 , 0.5 , 0.1 , 0.1,
#                    0.1, 0.1  , 0.4 , 0.4,
#                    0.1, 0.1  , 0.4 , 0.4),
#                  nrow=4, byrow=T)


# -> simple asymmetric transition matrix
#    should be clustered as {(1,2,3), (4), (5,6,7)}
#
# K <- matrix(c(0.5, 0.3, 0.2,   0,   0,   0,   0,
#               0.2, 0.6, 0.2,   0,   0,   0,   0,
#               0.3, 0.3, 0.4,   0,   0,   0,   0,
#               0  ,   0,0.05, 0.9,0.05,   0,   0,
#               0  ,   0,   0,   0, 0.6, 0.2, 0.2,
#               0  ,   0,   0,   0,0.35, 0.4,0.25,
#               0  ,   0,   0,   0, 0.3, 0.2, 0.5),
#             nrow=7, byrow=T)
# T <- rep(1,7)
# Tp <- rep(1,7)
#
# clusters <- asym.spectral.clustering(K, 3)
# wcut(K, T, Tp, clusters) 
