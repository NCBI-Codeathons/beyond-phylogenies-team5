
pairs_ti_clustrer <- function(input) {
#to convert pairwise data (e.g from TransPhylo output) into cluster data

#install.packages('igraph')
library(igraph)
#Addtifyverse
require (tidyverse)

#read the pairwise data file passed as argument
df_pairs <- read.table(file = input, sep = '\t', header = TRUE)
  
#df_pairs <- read.table(file = 'transphylo_demo_output.tsv', sep = '\t', header = TRUE)

df_pairs$connection <- 1

#transform into a matrix
m_pairs <- graph.data.frame(df_pairs, directed=FALSE)

#make adjacency matrix
matrix_pairs <- as_adjacency_matrix(m_pairs, names=TRUE, sparse=FALSE, attr="connection")
diag(matrix_pairs) <- 0

#write.csv(matrix_pairs,"matrix_pairs.csv")

#sorting the adjacency matrix into bottom-up clusters, in a way that optimizes the modularity score
community2 <- fastgreedy.community(as.undirected(graph.adjacency(matrix_pairs)))
sizes(community2)

#output dataframe with each sample's cluster ID written next to it
df_clusters <- as.data.frame(cbind(clusterID=community2$membership, sample=community2$names))
#
# make it wider format with number of members in cluster and details
df_clusters_rowwise <- df_clusters%>% group_by(clusterID) %>% summarize (N=n(),Clusters=paste (sample,collapse = ":"))
# retuns clusters two ways in a list
return(df_clusters.df_clusters_rowwise)
}
#write.csv(df_clusters,"cluster_assignments.csv", row.names=FALSE)
