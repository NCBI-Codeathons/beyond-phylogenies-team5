library(data.table)
library(tidyverse)
library(lubridate)
library(transcluster)

setwd("/Users/cc19/Documents/NCBI CODEATHON/DATA/SIMULATED_DATASET")

snp_dist<-as_tibble(fread("sars-cov-2_simulation_output.SNPDIST.tsv",sep="\t",header=TRUE)) %>%
  rename(strain1=`snp-dists 0.7.0`) 
snp_dist_final<-snp_dist %>% 
  gather(key="strain2",value="snp.dist",-strain1) %>% rowwise() %>%
  mutate(XX=ifelse(strain1>strain2, paste0(strain1,"/",strain2),paste0(strain2,"/",strain1) )) %>%
  group_by(XX) %>% dplyr::filter(row_number()==1)

sample_dates<-sapply( unique(c(snp_dist_final$strain1,snp_dist_final$strain2)),function(x){ str_split(x,"_")[[1]][[2]] } )
sample_data<-list2DF(list(sample.id=names(sample_dates),sample.dates=sample_dates)) %>%
  column_to_rownames("sample.id")

snp.matrix<-snp_dist %>% column_to_rownames("strain1")
dim(snp.matrix)

sample.ids<-colnames(snp.matrix)
sample.dates<-unname(sapply(sample_data[sample.ids,],function(x){  decimal_date(as.Date(x))  }))
##snp.matrix

sarscov2.mol.clock<-24/365

transcl.model <- createModel(sample.ids,sample.dates,snp.matrix)
transcl.model <- setParams(transcl.model,lambda=sarscov2.mol.clock)

#transcl.model <- setSNPThresholds(transcl.model, c(5,6,7))
#transcl.model <- setTransThresholds(transcl.model, c(7,8,9))

mySNPClusters <- makeSNPClusters(transcl.model,'model')
myTransClusters <- makeTransClusters(transcl.model,'model')

plotSNPClusters(mySNPClusters, vColor='orange', vSize=20, vFontSize=2, level=5, thick=3, labelOffset=0)
plotTransClusters(myTransClusters, vColor='orange', vSize=20, vFontSize=2, level=5, thick=3, labelOffset=0)

myModel <- setCutoffs(transcl.model)
plotTransClusters(myModel, vSize=4, vFontSize=0.5, level=7, thick=0.3, labelOffset=0)




