library(dplyr)
library(tidyr)
library(data.table) 
library(reshape2)

#read the data file with row-wise clusters (e.g. output from transcluster)
df_input <- t(read.table(file="transcluster_demo_output.csv", sep = ",", fill = T))
df_input <- df_input[,-1]

#assign column ID to each column
colnames(df_input) <- paste0("cluster#", 1:ncol(df_input))

df_clusters <- reshape2::melt(df_input)
df_clusters[df_clusters==""] <- NA
df_clusters <- na.omit(df_clusters)

write.csv(df_clusters,"cluster_assignments.csv", row.names=FALSE)
#output has all sample listed row-wise, along with their cluster ID