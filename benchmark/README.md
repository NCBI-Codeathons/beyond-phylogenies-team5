To benchmark and compare the output from three different tools (TransPhylo, outbreaker2, and transcluster). The idea is to get a common output format that can then be used to compare the three methods.

1. Function convert_transpylo (in pairs_to_cluster.R) takes TransPhylo output in pairwise form (e.g. transphylo_demo_output.tsv) and outputs cluster assignments alongside each sample.
2. Function convert_outbreaker (in outbreaker_to_clusters.R) takes outbreaker2 output in pairwise form (e.g. outbreaker2_demo_output.csv) as well as a probability threshold value (recommended: 0.72) and outputs cluster assignments alongside each sample.

