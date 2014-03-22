

# # # # # # #
# FUNCTIONS #
# # # # # # #

# The functions below are implementations of methods from Meila & Pentney
# asym.spectral.clustering - corresponds to the *BestWCut* algorithm
# wcut - corresponds to the quantity *WCut* defined in the paper
#

asym.spectral.clustering <- function(G, K){
  # following Meila & Pentney: Clustering by weighted cuts in directed graphs
  
  HB <- 0.5 * (2 * diag(nrow(G)) - G - t(G))
  X <- eigen(HB)$vectors[,tail(seq(nrow(HB)), K)]
  T <- t(apply(X, 1, function(x) x/sqrt(sum(x^2))))
  return(kmeans(T,K)$cluster)
}

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
