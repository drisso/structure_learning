# Code to reproduce the analyses of Nguyen et al. (2020)

In order to run the code in this repository, you need first to install the 
`zinbgraph` R package, available at https://github.com/drisso/zinbgraph.

## Simulation results

To reproduce simulation results, run `scripts/simulation_examples.R`.

Change the structure of the graph to simulate from, by changing the loaded file
in line 11. You can control the cardinality, number of simulations, and alpha
level by changing the values of lines 14-20.

To see example code of how to create the graph structure that we used for our
simulations, see `scripts/samplegraphs.R`.

## Real data

To reproduce results on real data:

1. run `scripts/learn_tf_real_data.R` to learn the structure of the graphs for different developmental stages.
2. run `scripts/geneset_enrichment.Rmd` to perform the Leiden community detection and gene set enrichment analysis.
3. run `scripts/TF_leiden_hiveplot.Rmd` to visualize the networks with hive plots.

Note that we included the count matrix in the `data/counts.rds` object for convenience.
The full, raw data are available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153730.
