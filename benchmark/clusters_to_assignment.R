require(tidyr, data.table, reshape2)

#read the clustered data file, transpose it, and drop the text about the threshold
df_input <- t(read.table(file="transcluster_demo_output.csv", sep = ",", fill = T))
df_input <- df_input[,-1]

#assign column ID to each column
colnames(df_input) <- paste0("cluster#", 1:ncol(df_input))

df_clusters <- reshape2::melt(df_input)