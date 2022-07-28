library(ape)
library(outbreaker2)
library(distcrete)
setwd('~/Projects/beyond-phylogenies-team5/outbreaker2/sim/')

#read in simulated data
seqs = read.FASTA('../../simulation/sim_123-23/sars-cov-2_simulation_output.fasta.tar.gz',type='DNA')

dates = strsplit(names(seqs),'_')
dates = as.Date(sapply(dates, "[[", 2))
minDate = min(dates)
dates = dates - minDate
names(dates) = names(seqs)
### generate serial interval distribution
# sim data, gamma distribution 
mu=5.2
sigma=1.72
shape = (mu/sigma)^2
scale = sigma^2/mu
si <- distcrete("gamma", shape = shape,scale = scale,interval = 1L, w = 0L)
si = si$d(1:15)
write.table(si,'sim_genint.csv')
pdf('sim_genint.pdf')
plot(si)
dev.off()

#also generate for MA real-world data
mu=2.9
sigma=1.9
shape = (mu/sigma)^2
scale = sigma^2/mu
si <- distcrete("gamma", shape = shape,scale = scale,interval = 1L, w = 0L)
si = si$d(1:15)
write.table(si,'realworld_genint.csv')

pdf('realworld_genint.pdf')
plot(si)
dev.off()


## assemble data for outbreaker
data <- outbreaker_data(dna = seqs, dates = dates, w_dens = si)

set.seed(1)
start_time <- Sys.time()
res <- outbreaker(data = data)
end_time <- Sys.time()
saveRDS(res, file = 'res_sim.rds')
runtime=end_time - start_time
print(runtime)


pdf('plot_posterior.pdf')
plot(res,"post")
dev.off()

pdf('plot_prior.pdf')
plot(res,"prior")
dev.off()

pdf('plot_rate.pdf')
plot(res,"mu")
dev.off()

pdf('mu_hist.pdf')
plot(res, "mu", "hist", burnin = 2000)
dev.off()

pdf('alpha_post_freq.pdf')
plot(res, type = "alpha", burnin = 2000)
dev.off()

library(htmltools)
save_html(plot(res, type = "network", burnin = 3000, min_support = 0.01),'test.html')

