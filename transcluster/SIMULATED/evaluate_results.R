source("../benchmark/pairs_to_cluster.R")
library(dplyr)
library(ggplot2)

true_pairs <- convert_transpylo("../simulation/sim_123-23/sim_123-23.pairs")
true_clusters <- n_distinct(true_pairs$clusterID)

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

output_files <- list.files("./", pattern="*.csv")
res <- lapply(output_files, function(x){
    read_input(x)
})
res_df <- do.call("rbind", res)
known_res_df <- res_df %>% filter(seq_id %in% true_pairs$sample)

nclusters_df <- known_res_df %>%
    group_by(threshold, threshold_val) %>%
    summarise(ncluster = n_distinct(cluster))


nclusters_df %>%
    ggplot(aes(threshold_val, ncluster)) + geom_line() + geom_point() + geom_hline(yintercept = true_clusters, linetype = "dashed", color="red") + facet_grid(threshold ~ .) +theme_bw() + scale_x_continuous(breaks=seq(1, 10))
ggsave("./simulated_results.pdf", w= 7.5, h = 10)
