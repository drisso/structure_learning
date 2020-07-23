library(slingshot)
library(clusterExperiment)

counts <- readRDS("extdata/finalTrajectory/counts_noResp_noMV.rds")
labels <- read.csv("extdata/finalTrajectory/Datta_annotation_HBC_lineage_newPCs.csv", row.names = 1)
labels <- labels[colnames(counts),]

load("extdata/regenK5_1_se_filtqc_idfiltyes_20190530_160100.Rda")

stopifnot(all(colnames(counts) %in% colnames(se_filtered)))

table(labels$tp, labels$cl)
table(labels$cl_names, labels$cl)

library(ggplot2)
ggplot(labels, aes(x=UMAP_1, y=UMAP_2, color=cl_names)) +
    geom_point()

labs <- readRDS("extdata/finalTrajectory/dattaCl_noResp_noMV.rds")
table(labs, labels$leiden_0_6_caleb20PCs)
labels$labs <- labs

ggplot(labels, aes(x=UMAP_1, y=UMAP_2, color=labs)) +
    geom_point()

stopifnot(all(rownames(labels) == colnames(counts)))

neuronal_counts <- assay(se_filtered)[,colnames(counts)]
stopifnot(all(colnames(neuronal_counts) == colnames(counts)))

saveRDS(neuronal_counts, file = "counts.rds")
