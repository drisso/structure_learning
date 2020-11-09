# Code to reproduce the analyses of Nguyen et al. (2020)

## Simulation results

To reproduce simulation results, run `scripts/simulation_examples.R`.

Change the structure of the graph to simulate from, by changing the loaded file
in line 11. You can control the cardinality, number of simulations, and alpha
level by changing the values of lines 14-20.

To see example code of how to create the graph structure that we used for our
simulations, see `scripts/samplegraphs.R`.

## Real data

To reproduce results on real data:

1. run `scripts/prepare_data.R` to preprocess the real dataset.
2. run `scripts/learn_tf_real_data.R` to learn the structure of the graphs for different developmental stages.
3. run `scripts/geneset_enrichment.Rmd` to perform the Leiden community detection and gene set enrichment analysis.
4. run `scripts/TF_leiden_hiveplot.Rmd` to visualize the networks with hive plots.