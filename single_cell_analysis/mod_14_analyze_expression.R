# clear existing variables
rm(list=ls())

library("tidyverse")
library("scales")
library("gridExtra")
library("hrbrthemes")
options(tibble.width = Inf)
library("Seurat")


set.seed(2)
t1 <- Sys.time()
rds_dir <- "rds"

# load data ===============

message("Loading data...")
t2 <- Sys.time()
x <- readRDS(file = file.path(rds_dir, "rds_13_data_filtered.rds"))
cat("\t"); print(Sys.time() - t2)


# Principal component analysis ===============

message("Applying PCA...")
n_pcs <- 40
############# modified code
x<- RunPCA(object      = x,
           features = VariableFeatures(object = x),
           do.print    = TRUE,
           pcs.compute = n_pcs, # store the first 40 principal components
           pcs.print   = 5,     # only print 5 of them
           genes.print = 5,     # print 5 genes for each PC
           use.imputed = FALSE,
           rev.pca     = FALSE)

# determine statistically significant PCs ===============

p1 <- ElbowPlot(x, ndims = 40)

print(p1)

# save data ===============

message("Saving data...")
t2 <- Sys.time()
saveRDS(x, file = "rds/rds_14_data_expression.rds")
cat("\t"); print(Sys.time() - t2)






message("Done.")
print(Sys.time() - t1)