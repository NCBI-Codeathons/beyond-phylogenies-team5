library(TransPhylo)
library(ape)

# note units of time here are days
# they will be converted to fraction of year when passed to sim function
Ne <- 10
start <- 2019.1
end <- 2019.2
# R0 = off.r*off.p/(1‐off.p)
off.r=4
pi <- 0.2
gen.mean <- 5.2
gen.sd <- 1.72
samp.mean <- 5.2
samp.sd <- 1.72
out<-"sim2"

#R <- off.r*off.p/(1.0-off.p)
#print(paste0("Using R value of ", R))

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

# run simulation (can take a while depending on params)
simu <- simulateOutbreak(neg=Ne*(gen.mean/365),
                         pi=pi,
                         off.r=off.r,
                         w.mean=gen.mean/365,
                         w.std=gen.sd/365,
                         ws.mean=samp.mean/365,
                         ws.std=samp.sd/365,
                         dateStartOutbreak=start,
                         dateT=end)
simu

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

# save simulation output as .rds
saveRDS(simu, file=paste0(out, ".rds"))

# output some stuff
print(paste0("Number of tips sampled: ", length(p$tip.label)))
print(paste0("Number of true pairs sampled:", nrow(pairs)))
print(paste0("Writing outputs with prefix: ", out))
print(paste0("Date last sample:"), dateLastSample(simu))
