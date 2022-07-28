library(outbreaker2)
library(ape)
library(distcrete)

# import sequence data
seq <- read.FASTA('', type = 'DNA')

# create vector for generation time distribution
mu = 2.9
sigma = 1.9
shape = (mu/sigma)^2
scale = sigma^2/mu
si = distcrete("gamma", shape = shape, scale = scale, interval = 1L, w = 0L)
si = si$d(1:15)

# get dates of case isolation by parsing gene sequence metadata
dates = strsplit(names(seq), '\\|')
dates = as.Date(sapply(dates, "[[", 3))
dates = dates - min(dates)
names(dates) = names(seq)

# process data and run outbreaker()
data <- outbreaker_data(dna = seq, dates = dates, w_dens = si)
set.seed(1)
start = Sys.time()
res <- outbreaker(data = data)
end = Sys.time()
cat('Time elapsed for outbreaker(): ', end - start)

# visualize results
pdf('sim_genint.pdf')
plot(si)
dev.off()

pdf('plot_posterior.pdf')
plot(res, "post")
dev.off()

pdf('plot_prior.pdf')
plot(res, "prior")
dev.off()

pdf('plot_rate.pdf')
plot(res, "mu")
dev.off()

pdf('mu_hist.pdf')
plot(res, "mu", "hist", burnin = 2000)
dev.off()

pdf('alpha_post_freq.pdf')
plot(res, type = "alpha", burnin = 2000)
dev.off()