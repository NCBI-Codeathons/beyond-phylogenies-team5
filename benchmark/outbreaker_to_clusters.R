#to convert pairwise transmission data (from outbreaker2 full_link_table) into cluster data

library(igraph)
library(plyr)
library(dplyr)

convert_outbreaker <- function(fname, thresh) {
  #read the pairwise data file with the associated posterior probabilities
  df_input <- read.table(file = fname)
  
  #subsetting by a threshold of posterior probabilities
  sub_df <- df_input[df_input$support >= thresh, ]
  
  #subsetting the columns with pairwise sample names and inferred generations
  cols_to_subset <- c("fromT","toT","generations")
  df_pairs <- df_o[,cols_to_subset]
  
  #transform into a matrix
  m_pairs <- graph.data.frame(df_pairs, directed=FALSE)
  
  #make adjacency matrix
  matrix_pairs <- as_adjacency_matrix(m_pairs, names=TRUE, sparse=FALSE, attr="generations")
  diag(matrix_pairs) <- 0
  
  #sorting the adjacency matrix into bottom-up clusters, in a way that optimizes the modularity score
  community2 <- fastgreedy.community(as.undirected(graph.adjacency(matrix_pairs)))
  sizes(community2)
  
  #output dataframe with each sample's cluster ID written next to it
  df_clusters <- as.data.frame(cbind(clusterID=community2$membership, sample=community2$names))
  
  return(df_clusters)
}

