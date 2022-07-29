options(rlib_downstream_check = FALSE)
suppressPackageStartupMessages({
  library(TransPhylo)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(khroma)
})

#################################################
# Functions
#################################################

get_tpr_fpr_pairs <- function(exp, obs, nseqs = 123){
    sort_combine <- function(x) {paste(sort(x), collapse=" ")}
    exp_pairs <- apply(exp, 1, sort_combine)
    obs_pairs <- apply(obs, 1, sort_combine)

    tp <- sum(obs_pairs %in% exp_pairs)
    fp <- sum(!(obs_pairs %in% exp_pairs))
    fn <- sum(!(exp_pairs %in% obs_pairs))
    tn <- (nseqs * (nseqs - 1) ) - (tp + fp + fn)

    tpr <- tp/(tp+fn)
    fpr <- fp/(fp + tn)

    yj <- (tp/(tp + fn)) + (tn/(tn + fp)) - 1
    res <- data.frame(
        tpr = tpr,
        fpr = fpr,
        yj = yj
    )
    return(res)
}

get_pairs <- function(wiw, thresh){
  to<-gsub(".*_", "", rownames(wiw)[which(wiw > thresh, arr.ind = TRUE)[, 1]])
  from<-gsub(".*_", "", colnames(wiw)[which(wiw > thresh, arr.ind = TRUE)[, 2]])
  pp<-wiw[wiw > thresh]
  pairs <- data.frame("from"=from, "to"=to, "prob"=pp)
  return(pairs)
}

# Trancluster
read_transcluster_input <- function(fname){
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

get_pairs_from_clusters <- function(clusters){
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

# Read in true pairs
true_pairs <- read.table("../simulation/sim_123-23/sim_123-23.pairs", header=T, sep="\t")

# TransPhylo: read outputs
mat <- readRDS("../simulation/sim_123-23/transphylo_output/transphylo_sim123-23_transphylo_matrix.rds")

# Get TPR and FPT for different clustering thresholds
results<-list()
for (thresh in seq(0, 1, 0.05)){
  pairs <- get_pairs(mat, thresh)
  pairs <- pairs[,1:2]
  stats <- get_tpr_fpr_pairs(true_pairs, pairs)
  stats["thresh"] <- thresh
  stats["tool"] <- "TransPhylo"
  results <- rbind(results, stats)
}

# Outbreaker2: read outputs
links <- read.table("../outbreaker2/sim/full_link_table.csv")
for (thresh in seq(0, 1, 0.05)){
    linksLoc <- links[links$support>=thresh,c('fromT','toT')]
    pairs <- pairs[,1:2]
    stats <- get_tpr_fpr_pairs(true_pairs, linksLoc)
    stats["thresh"] <- thresh
    stats["tool"] <- "Outbreaker2"
    results <- rbind(results, stats)
}

# Transcluster
output_files <- list.files("../transcluster/SIMULATED/", pattern="*.csv", full.names=T)

res <- lapply(output_files, function(x){
    read_transcluster_input(x)
})
res_df <- do.call("rbind", res)
res_df$seq <- gsub(".*_", "", res_df$seq)

for (group_df in split(res_df, list(res_df$threshold, res_df$threshold_val))){
  threshold_type<-group_df[1,"threshold"]
  thresh<-group_df[1,"threshold_val"]
  print(paste0(threshold_type," -- ", thresh))
  pairs <- get_pairs_from_clusters(group_df)
  stats <- get_tpr_fpr_pairs(true_pairs, pairs)
  stats["thresh"] <- thresh
  stats["tool"] <- paste("transcluster", threshold_type, collapse=": ")
  results <- rbind(results, stats)
}

# Plots
results %>%
    arrange(tpr, fpr) %>%
    ggplot(aes(fpr, tpr, color = tool)) + geom_abline(method="lm", linetype="dashed") + geom_line() + geom_point() + theme_bw() + scale_color_bright() + xlab("TPR") + ylab("FPR")
ggsave("./ROC.png", w= 7.5, h = 5)

results %>%
    arrange(tpr, fpr) %>%
    filter(!grepl("transcluster", tool)) %>%
    ggplot(aes(fpr, tpr, color = tool)) + geom_abline(method="lm", linetype="dashed") + geom_line() + geom_point() + theme_bw() + scale_color_bright() + xlab("TPR") + ylab("FPR")
ggsave("./ROC_subset.png", w= 7.5, h = 5)

results <- read.csv("./tpr_fpr_yj.csv")

# Plot Yj index
results %>%
    mutate(label = paste(tool, thresh, sep=": ")) %>%
    filter(thresh >= 0.5) %>%
    ggplot(aes(reorder(label, yj), yj, fill = tool)) + geom_col() + coord_flip() + theme_bw() + scale_fill_bright() + ylab("Youden Index") + xlab("Tool: threshold")
ggsave("./yj.png", w= 7.5, h = 10)
