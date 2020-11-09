library(igraph)
################# sample scalefree graph 
graph <- barabasi.game(n=10, m = 2, out.dist = NULL, out.seq = NULL, out.pref = FALSE, 
              directed=FALSE)
plot(graph)

SAdj <- as_adjacency_matrix(graph)
sum(SAdj)

# save(SAdj, file = "scalefreesample10.RData")

################# sample random graph 

g <-erdos.renyi.game(n=10, p=0.3, type=c("gnp"),
                     directed = FALSE, loops = FALSE)
SAdj <- as_adj(g)
plot(g)
sum(SAdj)/2


# save(SAdj, file = "randomsample10-03.RData")

################# sample hub graph 
library(XMRF)
SAdj <- XMRF.Sim(n = 100, p = 10, model = "LPGM", graph.type = "hub")$B
graph <- graph_from_adjacency_matrix(SAdj, mode =  "undirected", weighted = NULL,
                                     diag = TRUE, add.colnames = NULL, add.rownames = NA)
plot(graph)
sum(SAdj)

# save(SAdj, file = "hubsample100.RData")


###############3
