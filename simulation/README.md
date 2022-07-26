# Data simulation -- "known pairs" dataset

## Simulating transmission tree in TransPhylo
Execute script `sim_transphylo.R`:

```
Rscript --vanilla ./sim_transphylo.R
```

For now, parameters are hardcoded inside this file, but I can add argument parsing later on:
```
Ne <- 10
start <- 2019.1
end <- 2019.2
# R0 = off.r*off.p/(1‐off.p)
off.r=4
pi <- 0.1
gen.mean <- 5.2
gen.sd <- 1.72
samp.mean <- 5.2
samp.sd <- 1.72
out<-"sim"
```

These will be passed to the `TransPhylo::simulateOutbreak` function like so:
```
simu <- simulateOutbreak(neg=Ne*(gen.mean/365),
                         pi=pi,
                         off.r=off.r,
                         w.mean=gen.mean/365,
                         w.std=gen.sd/365,
                         ws.mean=samp.mean/365,
                         ws.std=samp.sd/365,
                         dateStartOutbreak=start,
                         dateT=end)
```

The script will generate several outputs:
- $out.pairs: True pairs formatted as .tsv
- $out.tre: Time tree from simulateOutbreak formatted as newick
- $out.rds: simulateOutbreak object saved as RDS


## Simulating COVID-like sequences using phastSim

Install [phastSim](https://github.com/NicolaDM/phastSim):
```
# dependencies
pip install numpy importlib-resources six ete3 biopython protobuf google

# install phastsim
pip install phastSim

# see github for alternate instructions (including conda):
# https://github.com/NicolaDM/phastSim
```

Simulate sequences using COVID reference and UNREST substitution model tuned for COVID-19. Note there is a lot more capability to tune the substitution model than we are using here.
```
phastSim --outpath ./ --reference MN908947.3.fasta --treeFile sim.tre --createPhylip --createFasta --seed 10
```

## Troubleshooting
