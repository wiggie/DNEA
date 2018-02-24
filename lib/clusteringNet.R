library(igraph)
set.seed(1)
getConsensusMatrix <- function(cl){
  K <- length(cl)
  p <- length(cl[[1]]$membership)
  D <- matrix(0, p, p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      for (k in 1:K){
        tmp = cl[[k]]$membership
        D[i,j] <- D[i,j] + (length(unique(tmp[c(i,j)]))==1)
      }
    }
  }
  D <- D + t(D)
  D <- D/K
  diag(D) <- rep(1, p)
  return(D)
}

## ensemble corresponds to using all clustering techniques to find a consensus clustering
## Input: 
## graph: an igraph object
## K: the number of clusters
## tau: the threshold when evaluating the consensus matrix
run_consensus_cluster <- function(graph, K=10, tau=0.5, 
                                  method=c("infomap","lpm","walktrap","ensemble"), 
                                  maxIter=5){
  method <- match.arg(method)
  if (method=="ensemble"){
    cl <- list()
    set.seed(1)
    cl[[1]] <- cluster_edge_betweenness(graph, weights = NULL)
    cl[[2]] <- cluster_fast_greedy(graph, weights = NULL)
    cl[[3]] <- cluster_infomap(graph, e.weights = NULL)
    cl[[4]] <- cluster_label_prop(graph, weights = NULL)
    cl[[5]] <- cluster_leading_eigen(graph, weights = NULL)
    cl[[6]] <- cluster_louvain(graph, weights = NULL)
    cl[[7]] <- cluster_walktrap(graph, weights = NULL)
  } else {
    cl <- vector("list", K)
    for (k in 1:K){
      set.seed(k)
      if (method=="lpm"){
        cl[[k]] <- cluster_label_prop(graph, weights = NULL)
      } else if (method=="infomap"){
        cl[[k]] <- cluster_infomap(graph, e.weights = NULL)
      } else if (method=="walktrap"){
        cl[[k]] <- cluster_walktrap(graph, weights = NULL, steps=k)
      } 
    }
  }
  D <- getConsensusMatrix(cl)
  iter <- 0
  while(length(table(D))>2 && iter<maxIter){
    diag(D) <- 0
    thresholded_D <- D*(D>tau)
    
    Dgraph <- graph.adjacency(thresholded_D, mode="undirected", weighted = TRUE)
    if (method=="ensemble"){
      dcl <- list()
      set.seed(iter)
      dcl[[1]] <- cluster_edge_betweenness(Dgraph, weights = E(Dgraph)$weight)
      dcl[[2]] <- cluster_fast_greedy(Dgraph, weights = E(Dgraph)$weight)
      dcl[[3]] <- cluster_infomap(Dgraph, e.weights = E(Dgraph)$weight)
      dcl[[4]] <- cluster_label_prop(Dgraph, weights = E(Dgraph)$weight)
      dcl[[5]] <- cluster_leading_eigen(Dgraph, weights = E(Dgraph)$weight)
      dcl[[6]] <- cluster_louvain(Dgraph, weights = E(Dgraph)$weight)
      dcl[[7]] <- cluster_walktrap(Dgraph, weights = E(Dgraph)$weight)
    } else {
      dcl <- vector("list",K)
      for (k in 1:K){
        set.seed(k)
        if (method=="lpm"){
          dcl[[k]] <- cluster_label_prop(Dgraph, weights = E(Dgraph)$weight)
        } else if (method=="infomap"){
          dcl[[k]] <- cluster_infomap(Dgraph, e.weights = E(Dgraph)$weight)
        } else if (method=="walktrap"){
          dcl[[k]] <- cluster_walktrap(Dgraph, weights = NULL, steps=k)
        }
      }
    }
    D <- getConsensusMatrix(dcl)
    iter <- iter + 1
  }
  new.order <- order(dcl[[1]]$membership)
  return(list(dcl=dcl[[1]]$membership,D=D,order=new.order,iter=iter))
}
