library(TransPhylo)
library(ape)
library(lattice)
library(coda)

set.seed(0)

simu <- readRDS("~/beyond-phylogenies-team5/simulation/sim.rds")

phy<-read.tree(file="~/beyond-phylogenies-team5/simulation/sim.tre")

ptree<-ptreeFromPhylo(phy,dateLastSample=dateLastSample(simu))


w.mean <- 5
w.std <- 2
dateT=Inf
mcmc=100000
thin=100
out="transphylo_test1"

res<-inferTTree(ptree,
                w.mean=w.mean,
                w.std=w.std,
                dateT=dateT,
                mcmcIterations=mcmc,
                thinning=thin,
                optiStart=1, 
                updateNeg = TRUE,
                updateOff.r = TRUE,
                updatePi = TRUE,
                updateOff.p = TRUE)

# plot mcmc 
pdf(paste0(out,"_mcmcout.pdf"))
plot(res)
dev.off()

# calculate ess
mcmc<-convertToCoda(res)
ess<-effectiveSize(mcmc)
print(ess)

#get tree
med=medTTree(res)
ttree=extractTTree(med)

# pairwise transmission matrix 
wiw <- computeMatWIW(res)
pdf(paste0(out,"_levelplot.pdf"))
lattice::levelplot(wiw,xlab='',ylab='')
dev.off()

# write wiw
to<-rownames(wiw)[which(wiw > 0.95, arr.ind = TRUE)[, 1]]
from<-colnames(wiw)[which(wiw > 0.95, arr.ind = TRUE)[, 2]]
pp<-which(wiw > 0.95, arr.ind = TRUE)
pairs <- data.frame("from"=from, "to"=to)
write.table(pairs,
            file=paste0(out, ".pairs"),
            col.names=TRUE,
            row.names=FALSE,
            sep="\t",
            quote=FALSE)

# pairwise # of intermediates 
