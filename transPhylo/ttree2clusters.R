library(igraph)
library(dbscan)
library(TransPhylo)
library(dplyr)

tp.out <- "~/beyond-phylogenies-team5/simulation/sim_20-2/sim_20-2.rds"
int.host.min <- 0
int.host.max <- 5
burn=0.2
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
  dist<-shortest.paths(g, v=ttree$nam, to=ttree$nam)
  dist[dist!=0] <- dist[dist!=0]-1 
}

res<-readRDS(tp.out)
ttree<-extractTTree(res)
g<-graph_from_ttree(ttree)

pdf(out)
plot_transmission_graph(g)
dev.off()


