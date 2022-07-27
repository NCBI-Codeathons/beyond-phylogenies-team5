library(plyr)

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
    df["threshold_val"] <- snp_threshold[2,1]
    return(df)
}

output_files <- list.files("./", pattern="*.csv")
res <- lapply(output_files, function(x){
    read_input(x)
})
res_df <- do.call("rbind", res)

true_pairs <- read.table("../simulation/sim_123-23/sim_123-23.pairs", header=T, sep="\t")
