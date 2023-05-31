# clear existing variables
rm(list=ls())

library("tidyverse")
library("scales")
library("gridExtra")
library("hrbrthemes")
options(tibble.width = Inf)
library("Seurat")


set.seed(2)
rds_dir <- "rds"
t1 <- Sys.time()

cell_to_sample_name <- function(x) {
  sub("(.*)_(.*)", "\\2", x)
}

# load data ===============
message("Loading data...")
x <- readRDS(file = "rds/rds_14_data_expression.rds")

temp <- readRDS(file = "rds/rds_14_data_expression.rds")

# find clusters ===============

message("Clustering data...")
# save.SNN=T saves the SNN so that the  SLM algorithm can be rerun
# using the same graph, but with a different resolution value (see docs
# for full details)
t2 <- Sys.time()
dims_use <- 1:20
# first compute KNN graph
x <- FindNeighbors(
  object = x,
  dims = 1:20, reduction = "pca"   # specify the dimensions to use here
)


x <- FindClusters(object = x, resolution = 0.35)
# this identified 24 clusters

# dimension reduction ===============
x <- RunTSNE(x,
             dims.use = dims_use,
             do.fast = TRUE)
TSNEPlot(x)


# save data ===============

message("Saving data...")
t2 <- Sys.time()
saveRDS(x, file = "rds/rds_15_data_clusters.rds")
cat("\t"); print(Sys.time() - t2)






message("Done.")
print(Sys.time() - t1)