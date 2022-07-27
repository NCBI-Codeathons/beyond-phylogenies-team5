library(TransPhylo)
library(ape)
library(lattice)
library(coda)

set.seed(12322)

setwd("~/beyond-phylogenies-team5/simulation/sim_20-2/")

simu <- readRDS("sim_20-2.rds")

phy<-read.tree(file="sim_20-2.tre")

ptree<-ptreeFromPhylo(phy,dateLastSample=dateLastSample(simu))


w.mean <- 5.2
w.std <- 1.72
dateT=2019.2
mcmc=500000
thin=1000
burn=0.2
out="~/beyond-phylogenies-team5/transPhylo/transphylo_sim20-2"

res<-inferTTree(ptree,
                w.mean=w.mean/365,
                w.std=w.std/365,
                dateT=dateT,
                mcmcIterations=mcmc,
                thinning=thin,
                optiStart=2, 
                updateNeg = TRUE,
                updateOff.r = TRUE,
                updatePi = TRUE,
                updateOff.p = FALSE)

# plot mcmc 
pdf(paste0(out,"_mcmcout.pdf"))
plot(res)
dev.off()

# calculate ess
mcmc<-convertToCoda(res)
ess<-effectiveSize(mcmc)
print(ess)

#get tree
med=medTTree(res, burnin=burn)
ttree=extractTTree(med)

cons.ttree <- consTTree(res, burnin=burn)

# pairwise transmission matrix 
wiw <- computeMatWIW(res, burnin=burn)
pdf(paste0(out,"_levelplot.pdf"))
lattice::levelplot(wiw,xlab='',ylab='')
dev.off()

# write wiw
to<-rownames(wiw)[which(wiw > 0.95, arr.ind = TRUE)[, 1]]
from<-colnames(wiw)[which(wiw > 0.95, arr.ind = TRUE)[, 2]]
pp<-wiw[wiw > 0.95]
pairs <- data.frame("from"=from, "to"=to, "prob"=pp)
write.table(pairs,
            file=paste0(out, ".pairs"),
            col.names=TRUE,
            row.names=FALSE,
            sep="\t",
            quote=FALSE)

# pairwise # of intermediates

saveRDS(res, paste0(out, "_transphylo_output.rds"))
saveRDS(wiw, paste0(out, "_transphylo_matrix.rds"))
