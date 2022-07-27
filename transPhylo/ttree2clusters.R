library(igraph)
library(TransPhylo)
library(dplyr)

tp.out <- "~/beyond-phylogenies-team5/simulation/sim_20-2/transphylo_sim20-2_transphylo_output.rds"
int.host.min <- 0
int.host.max <- 5
burn=0.2
threshold=2
out <- "~/beyond-phylogenies-team5/simulation/sim_20-2/transPhylo_sim_20-2.graph.pdf"

graph_from_ttree <- function(ttree){
  tips <- ttree$nam
  t <- data.frame(ttree$ttree)
  colnames(t)<- c("inf.date", "samp.date", "from")
  t<-t %>% mutate("to" = row_number()) 
  t<-t %>% dplyr::select(c("from", "to"))
  t$nodeType<-"intermediate"
  t$name<-t$to
  t[1:length(tips),]$nodeType <- "observed"
  vertices <- t %>% select(name, nodeType)
  edges <- t %>% select(from, to) %>% filter(from!=0)
  g<-graph_from_data_frame(edges, directed=TRUE, vertices=vertices)
  return(g)
}

plot_transmission_graph <- function(g){
  plot(g, edge.arrow.size=.5,
       vertex.color=c( "red", "black")[1+(V(g)$nodeType=="intermediate")],
       vertex.size=c(5, 1)[1+(V(g)$nodeType=="intermediate")], 
       vertex.label.dist=1)
}

cluster_transmission_graph <- function(g, threshold){
  # get dist mat as # theoretical hosts between pairs
  tips <- V(g)[V(g)$nodeType!="intermediate"]$name
  dist<-shortest.paths(g, v=tips, to=tips)
  dist[dist!=0] <- dist[dist!=0]-1 
  adjacency <- dist
  adjacency[dist<=threshold] <- 1
  adjacency[dist>threshold] <- 0
  diag(adjacency) <- 0
  plot(as.undirected(graph.adjacency(adjacency)))
  clusters<-fastgreedy.community(as.undirected(graph.adjacency(adjacency)))
  return(clusters)
}

res<-readRDS(tp.out)
med<-medTTree(res, burnin=burn)
ttree<-extractTTree(med)
g<-graph_from_ttree(ttree)

pdf(out)
plot_transmission_graph(g)
dev.off()

clust<-cluster_transmission_graph(g, 2)

df_clusters <- as.data.frame(cbind(clusterID=clust$membership, sample=clust$names))

df_clusters_rowwise <- df_clusters%>% group_by(clusterID) %>% summarize (N=n(),Clusters=paste (sample,collapse = ":"))

write.csv(df_clusters,
          paste0(out,"_cluster_assignments.csv"), 
          row.names=FALSE,
          quote=FALSE)

write.csv(df_clusters_rowwise,
          paste0(out,"_clusters.csv"), 
          row.names=FALSE,
          quote=FALSE)
