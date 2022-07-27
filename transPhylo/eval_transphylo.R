library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

# get pairs above some posterior probability threshold 
get_pairs <- function(wiw, thresh){
  to<-gsub(".*_", "", rownames(wiw)[which(wiw > thresh, arr.ind = TRUE)[, 1]])
  from<-gsub(".*_", "", colnames(wiw)[which(wiw > thresh, arr.ind = TRUE)[, 2]])
  pp<-wiw[wiw > thresh]
  pairs <- data.frame("from"=from, "to"=to, "prob"=pp)
  return(pairs)
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

# need to make a version for clusters 
#get_stats_clusters <- function(exp, obs\){
#  
#}

true_pairs <- read.table("~/beyond-phylogenies-team5/simulation/sim_123-23/sim_123-23.pairs",
                         header=T, 
                         sep="\t")

mat <- readRDS("~/beyond-phylogenies-team5/simulation/sim_123-23/transphylo_sim123-23_transphylo_matrix.rds")

results<-list()
for (thresh in seq(0.4, 1.0, 0.02)){
  pairs <- get_pairs(mat, thresh)
  stats <- get_stats_pairs(true_pairs, pairs)
  stats["prob.thresh"] <- thresh
  results <- rbind(results, stats)
}
results <- rbind.fill(results)

# make plots 
s <- results %>% select(tpr, fdr, fnr, precision, prob.thresh)
sm <- melt(s, id="prob.thresh")
pdf("~/beyond-phylogenies-team5/simulation/sim_123-23/transphylo_sim123-23_eval-plot.pdf")
ggplot(sm, aes(x=prob.thresh, y=value, color=variable)) +
  theme_minimal() + 
  geom_line(lwd=2)
dev.off()
  


