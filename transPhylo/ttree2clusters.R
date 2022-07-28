#!/usr/local/bin/Rscript

suppressPackageStartupMessages({
  library(igraph)
  library(TransPhylo)
  library(dplyr)
  library(GetoptLong)
})

#################################################
# Command-line arg parsing
#################################################

burn=0.2
threshold=2
out <- "out"

GetoptLong(
  "rds=s", "TransPhylo output (.rds output from run_transphylo.R)",
  "out=s", "Output file prefix",
  "burn=f", "Proportion of MCMC iterations to discard as burn-in",
  "thresh=f", "Number of intermediate hosts to build clusters"
)

#################################################
# Functions
#################################################

graph_from_ttree <- function(ttree){
  tips <- gsub(".*_", "", ttree$nam)
  t <- data.frame(ttree$ttree)
  colnames(t)<- c("inf.date", "samp.date", "from")
  t<-t %>% mutate("to" = row_number()) 
  t<-t %>% dplyr::select(c("from", "to"))
  t$nodeType<-"intermediate"
  t$from_name<-paste0("_",t$from)
  t$to_name<-paste0("_", t$to)
  for (i in 1:length(tips)){
    name<-tips[[i]]
    t[t$from==as.integer(i),"from_name"] <- name
    t[t$to==as.integer(i),"to_name"] <- name
  }
  t[1:length(tips),]$nodeType <- "observed"
  vertices <- t %>% select(to_name, nodeType)
  edges <- t %>% select(from_name, to_name) %>% filter(from_name!="_0")
  g<-graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
  return(g)
}

plot_transmission_graph <- function(g){
  plot(g, edge.arrow.size=.5,
       vertex.color=c( "red", "black")[1+(V(g)$nodeType=="intermediate")],
       vertex.size=c(5, 1)[1+(V(g)$nodeType=="intermediate")], 
       vertex.label.dist=1)
}

cluster_transmission_graph <- function(g, threshold, plot=FALSE, out="clusters"){
  # get dist mat as # theoretical hosts between pairs
  tips <- V(g)[V(g)$nodeType!="intermediate"]$name
  dist<-shortest.paths(g, v=tips, to=tips)
  dist[dist!=0] <- dist[dist!=0]-1 
  adjacency <- dist
  adjacency[dist<=threshold] <- 1
  adjacency[dist>threshold] <- 0
  diag(adjacency) <- 0
  if (plot == TRUE){
    pdf(paste0(out, "_clusterGraph.pdf"))
    plot(as.undirected(graph.adjacency(adjacency)))
    dev.off()
  }
  clusters<-fastgreedy.community(as.undirected(graph.adjacency(adjacency)))
  return(clusters)
}

#################################################
# Main
#################################################

print("Reading inputs...")
res<-readRDS(rds)

print("Extracting transmission tree...")
med<-medTTree(res, burnin=burn)
ttree<-extractTTree(med)

print("Building transmission graph...")
g<-graph_from_ttree(ttree)

pdf(paste0(out, "_tGraph.pdf"))
plot_transmission_graph(g)
dev.off()

print("Identifying clusters in transmission graph...")
clust<-cluster_transmission_graph(g, 2)

print("Writing outputs as _cluster_assignments.csv and _clusters.csv...")
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

print("Done!")
