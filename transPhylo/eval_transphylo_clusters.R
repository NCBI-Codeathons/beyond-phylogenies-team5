#!/usr/local/bin/Rscript

suppressPackageStartupMessages({
  library(igraph)
  library(TransPhylo)
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(GetoptLong)
})

#################################################
# Command-line arg parsing
#################################################

burn=0.2
min<-1
max<-10
out <- "out"

GetoptLong(
  "rds=s", "Path to TransPhylo output (.rds output from run_transphylo.R)",
  "truePairs=s", "Path to known .pairs file",
  "out=s", "Output file prefix",
  "burn=f", "Proportion of MCMC iterations to discard as burn-in",
  "min=f", "Minimum # of intermediate hosts threshold to test",
  "max=f", "Maximum # of intermediate hosts threshold to test"
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

# returns some stats 
# tp = # of true positives 
# fn = # of false negatives (real positives not inferred)
# fp = # of false positives (inferred pos not in true pairs)
# tpr = Sensitivity or true positive rate = TP / (TP+FN)
# fdr = false discover rate = FP / (TP+FP)
# fnr = miss rate or false negative rate = FN / (TP+FN)
# precision = precision or positive predictive value = TP / (TP+FP)
# tp.prop = TP / total expected positives
# fn.prop = FN / total expected positives
# fp.prop = FP / total inferred positives
# exp.pos = expected positives 
# obs.pos = observed positives
get_stats_pairs <- function(exp, obs, bidirectional=TRUE){
  fn<-0
  tp<-0
  fp<-0
  # get tp and fn
  for(i in 1:nrow(exp)){
    if (exp[i,1] %in% obs[,1]) {
      if (exp[i,2] %in% obs[obs[,1]==exp[i,1],2]){
        tp <- tp+1
        next
      }
    }else{
      if (bidirectional==TRUE){
        if (exp[i,1] %in% obs[,2]) {
          if (exp[i,2] %in% obs[obs[,2]==exp[i,1],1]){
            tp <- tp+1
            next
          }
        }
      }
    }
  }
  fp <- nrow(obs)-tp
  fn <- nrow(exp)-tp
  fnr <- fn/nrow(exp)
  sens <- tp / (tp+fn)
  fdr <- fp / (tp+fp)
  prec <- tp / (tp+fp)
  ret = data.frame(
    "fp" = fp,
    "tp" = tp,
    "fn" = fn,
    "exp.pos" = nrow(exp),
    "obs.pos" = nrow(obs),
    "tpr" = sens,
    "fdr" = fdr,
    "fnr" = fnr,
    "precision" = prec
  )
  return(ret)
}

get_pairs_from_clusters <- function(clusters){
  df_clusters <- as.data.frame(cbind(clusterID=clust$membership, sample=clust$names))
  df_clusters_rowwise <- df_clusters%>% group_by(clusterID) %>% summarize (N=n(),Clusters=paste (sample,collapse = ":"))
  from<-list()
  to<-list()
  for (i in 1:nrow(df_clusters_rowwise)){
    if (as.integer(df_clusters_rowwise[i,"N"]) > 1){
      inds <- strsplit(as.character(df_clusters_rowwise[i,"Clusters"]), ":")[[1]]
      pairs=combn(inds,2)
      for (j in 1:ncol(pairs)){
        from <- c(from, as.integer(pairs[1,j]))
        to <- c(to, as.integer(pairs[2,j]))
      }
    }
  }
  pairs <- data.frame("from"=unlist(from), "to"=unlist(to))
  return(pairs)
}


#################################################
# Main
#################################################

# read outputs and build graph
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

# evaluate for diff thresholds 
print("Reading true pairs...")
true_pairs <- read.table(truePairs,
                         header=T, 
                         sep="\t")


# try different clustering thresholds and get results 
print("Testing across specified clustering thresholds...")
results<-list()
for (thresh in seq(min, max, 1)){
  print(thresh)
  clust<-cluster_transmission_graph(g, thresh, plot=FALSE)
  pairs <- get_pairs_from_clusters(clust)
  if(nrow(pairs) < 1){
    print("No observed pairs... Skipping")
    next
  }
  stats <- get_stats_pairs(true_pairs, pairs)
  stats["thresh"] <- thresh
  results <- rbind(results, stats)
}
results <- rbind.fill(results)


# make plots 
print("Outputting results as _clustEval.pdf and _clustEval.tsv")
s <- results %>% select(tpr, fdr, fnr, precision, thresh)
sm <- melt(s, id="thresh")
pdf(paste0(out, "_clustEval.pdf"))
ggplot(sm, aes(x=thresh, y=value, color=variable)) +
  theme_minimal() + 
  geom_line(lwd=2)
dev.off()

write.table(results,
            paste0(out, "_clustEval.tsv"),
            col.names=TRUE,
            row.names=FALSE,
            sep="\t",
            quote=FALSE)

print("Done!")

