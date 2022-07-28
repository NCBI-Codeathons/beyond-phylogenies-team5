#!/usr/local/bin/Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(GetoptLong)
})

#################################################
# Command-line arg parsing
#################################################

out <- "out"

GetoptLong(
  "dir=s", "Directory containing transcluster outputs",
  "truePairs=s", "Path to known .pairs file",
  "out=s", "Output file prefix"
)

#################################################
# Functions
#################################################

read_input <- function(fname){
  lines <- readLines(fname)
  ctr <- 0
  rows <- lapply(lines, function(x){
    seq_ids <- strsplit(x, ",")[[1]]
    tmp <- data.frame(seq = seq_ids)
    tmp$cluster <- ctr 
    ctr <<- ctr + 1
    tmp
  })
  
  df <- do.call("rbind", rows)
  snp_threshold <- df[1:2,]
  df <- df[3:nrow(df),]
  df["threshold"] <- snp_threshold[1,1]
  df["threshold_val"] <- as.numeric(snp_threshold[2,1])
  df$seq_id <- df$seq %>%
    strsplit( "_" ) %>%
    sapply( "[", 3 )
  return(df)
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

get_pairs_from_clusterMembership <- function(clusters){
  from<-list()
  to<-list()
  for (clust in split(clusters, clusters$cluster)){
    if (nrow(clust) < 2){
      next
    }
    pairs=combn(clust$seq,2)
    for (j in 1:ncol(pairs)){
      from <- c(from, as.integer(pairs[1,j]))
      to <- c(to, as.integer(pairs[2,j]))
    }
  }
  pairs <- data.frame("from"=unlist(from), "to"=unlist(to))
  return(pairs)
}


#################################################
# Main
#################################################

print("Reading known pairs...")
true_pairs <- read.table(truePairs,
                         header=T, 
                         sep="\t")

print(paste0("Reading all transcluster outputs in", dir, "..."))
output_files <- list.files(dir, pattern="*.csv", full.names=T)

res <- lapply(output_files, function(x){
    read_input(x)
})
res_df <- do.call("rbind", res)
res_df$seq <- gsub(".*_", "", res_df$seq)

print("Testing across supplied clustering thresholds...")
results<-list()
for (group_df in split(res_df, list(res_df$threshold, res_df$threshold_val))){
  threshold_type<-group_df[1,"threshold"]
  thresh<-group_df[1,"threshold_val"]
  print(paste0(threshold_type," -- ", thresh))
  pairs <- get_pairs_from_clusterMembership(group_df)
  if(nrow(pairs) < 1){
    print("No observed pairs... Skipping")
    next
  }
  stats <- get_stats_pairs(true_pairs, pairs)
  stats["threshold"] <- thresh
  stats["threshold_type"] <- threshold_type
  results <- rbind(results, stats)
}
results <- rbind.fill(results)

nclusters_df <- known_res_df %>%
    group_by(threshold, threshold_val) %>%
    summarise(ncluster = n_distinct(cluster))


nclusters_df %>%
    ggplot(aes(threshold_val, ncluster)) + geom_line() + geom_point() + geom_hline(yintercept = true_clusters, linetype = "dashed", color="red") + facet_grid(threshold ~ .) +theme_bw() + scale_x_continuous(breaks=seq(1, 10))
ggsave("./simulated_results.pdf", w= 7.5, h = 10)
