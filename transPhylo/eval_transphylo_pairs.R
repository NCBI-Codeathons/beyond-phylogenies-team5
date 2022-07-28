#!/usr/local/bin/Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(reshape2)
  library(GetoptLong)
})

#################################################
# Command-line arg parsing
#################################################

minP<-0.5
maxP<-1.0
step<-0.05
out<-"out"

GetoptLong(
  "mat=s", "Path to pairwise transmission prob matrix (.rds output from run_transphylo.R)",
  "truePairs=s", "Path to known .pairs file",
  "out=s", "Output file prefix",
  "minP=f", "Minimum posterior probability threshold to test",
  "maxP=f", "Maximum posterior probability threshold to test",
  "step=f", "Step size for posterior probability range"
)

#################################################
# Functions
#################################################

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

#################################################
# Main
#################################################
print("Reading inputs...")
true_pairs <- read.table(truePairs,
                         header=T, 
                         sep="\t")

mat <- readRDS(mat)

print("Testing across specified clustering thresholds...")
results<-list()
for (thresh in seq(minP, maxP, step)){
  print(thresh)
  pairs <- get_pairs(mat, thresh)
  if(nrow(pairs) < 1){
    print("No observed pairs... Skipping")
    next
  }
  stats <- get_stats_pairs(true_pairs, pairs)
  stats["prob.thresh"] <- thresh
  results <- rbind(results, stats)
}
results <- rbind.fill(results)

# make plots 
print("Outputting results as _clustEval.pdf and _clustEval.tsv")

s <- results %>% select(tpr, fdr, fnr, precision, prob.thresh)
sm <- melt(s, id="prob.thresh")
pdf(paste0(out, "_pairEval.pdf"))
ggplot(sm, aes(x=prob.thresh, y=value, color=variable)) +
  theme_minimal() + 
  geom_line(lwd=2)
dev.off()
  
write.table(results,
            paste0(out, "_pairEval.tsv"),
            col.names=TRUE,
            row.names=FALSE,
            sep="\t",
            quote=FALSE)

print("Done!")
