library(ape)
library(outbreaker2)
setwd('~/Projects/beyond-phylogenies-team5/outbreaker2/')

col <- "#6666cc"

## initially we'll use the outbreaker2 test data
#sequence data stored in binary format
fake_outbreak

pdf('gen_time_dist.pdf')
### plot of generation time/serial interval distribution
plot(fake_outbreak$w, type = "h", xlim = c(0, 6), 
     lwd = 30, col = col, lend = 2, 
     xlab = "Days after infection", 
     ylab = "p(new case)", 
     main = "Generation time distribution")
dev.off()

# now let's actually do some analysis

dna <- fake_outbreak$dna
dates <- fake_outbreak$sample
ctd <- fake_outbreak$ctd
w <- fake_outbreak$w
data <- outbreaker_data(dna = dna, dates = dates, ctd = ctd, w_dens = w)

set.seed(1)
start_time <- Sys.time()
res <- outbreaker(data = data)
end_time <- Sys.time()
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

### now trying some more custom parameters 
config2 <- create_config(n_iter = 3e4,
                         sample_every = 20,
                         init_tree ="star",
                         move_kappa = FALSE,
                         prior_mu = 10)
set.seed(1)
start_time <- Sys.time()
res2 <- outbreaker(data, config2)
end_time <- Sys.time()
runtime=end_time - start_time
print(runtime)

library(htmltools)
save_html(plot(res2, type = "network", burnin = 3000, min_support = 0.01),'test.html')



