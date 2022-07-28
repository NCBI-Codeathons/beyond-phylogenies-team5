#!/usr/local/bin/Rscript

suppressPackageStartupMessages({
  library(TransPhylo)
  library(ape)
  library(lubridate)
  library(GetoptLong)
})

#################################################
# Command-line arg parsing
#################################################

ne <- 10
start <- 2019.1
end <- 2019.2
r <- 4
pi <- 0.1
genMean <- 5.2
genSD <- 1.72
sampMean <- 5.2
sampSD <- 1.72
nSampled=200
out<-"sim"
seed=321321

GetoptLong(
  "ne=i", "Intrahost effective population size (Ne)",
  "start=f", "Simulation start date",
  "end=f", "Simulation end date",
  "r=i", "R number",
  "pi=f", "Probability of an infected being sampled",
  "genMean=f", "Generation time distribution mean (days)",
  "genSD=f", "Generation time distribution stdev (days)",
  "sampMean=f", "Sampling time distribution mean (days)",
  "sampSD=f", "Sampling time distribution stdev (days)",
  "nSampled=i", "Number of individuals sampled",
  "out=s", "Output file prefix",
  "seed=i", "Random number seed"
)
  
#################################################
# Functions 
#################################################

# function to extract true sampled transmission pairs
extractTransmissionPairs <- function(ttree) {
  tips <- ttree$nam
  t <- data.frame(ttree$ttree)
  colnames(t)<- c("inf.date", "samp.date", "anc")
  from<-c()
  to<-c()
  for (tip in tips){
    a<-t[as.integer(tip),"anc"]
    if (as.character(a) %in% tips){
      from<-c(from, a)
      to <-c(to, tip)
    }
  }
  ret<-data.frame("from"=from, "to"=to)
  return(ret)
}

#function to annotate tip labels from ttree table
extractTipLabels <- function(ttree, tips){
  t <- data.frame(ttree$ttree)
  colnames(t)<- c("inf.date", "samp.date", "anc")
  ret<-c()
  for (tip in tips){
    index<-as.integer(tip)
    if (! is.na(t$samp.date[index])){
      d1<-format(as.Date(date_decimal(t[index,1])))
      d2<-format(as.Date(date_decimal(t[index,2])))
      nn <- paste0(d1, "_", d2, "_", tip)
      ret<-c(ret, nn)
    }
  }
  return(ret)
}

#################################################
# Main
#################################################
set.seed(seed)

print("Starting simulation...")
# run simulation (can take a while depending on params)
simu <- simulateOutbreak(neg=ne*(genMean/365),
                         pi=pi,
                         off.r=r,
                         w.mean=genMean/365,
                         w.std=genSD/365,
                         ws.mean=sampMean/365,
                         ws.std=sampSD/365,
                         dateStartOutbreak=start,
                         dateT=end,
                         nSampled=nSampled)
print("Finished simulation!")
simu

print("Extracting true pairs and writing outputs...")
ttree<-extractTTree(simu)
ptree<-extractPTree(simu)

# extract and write true pairs
pairs<-extractTransmissionPairs(ttree)
write.table(pairs,
            file=paste0(out, ".pairs"),
            col.names=TRUE,
            row.names=FALSE,
            sep="\t",
            quote=FALSE)

# format tips as "infected.date_sampled.date_index"
p<-phyloFromPTree(ptree)
p$tip.label <- extractTipLabels(ttree, p$tip.label)
write.tree(p, paste0(out, ".tre"))

saveRDS(ttree, paste0(out, ".ttree.rds"))

# save simulation output as .rds
saveRDS(simu, file=paste0(out, ".rds"))

# output some stuff
print(paste0("Number of tips sampled: ", length(p$tip.label)))
print(paste0("Number of true pairs sampled:", nrow(pairs)))
print(paste0("Writing outputs with prefix: ", out))

print("Done!")
