# clear existing variables
rm(list=ls())

library("tidyverse")
library("scales")
library("gridExtra")
#install.packages("hrbrthemes")
library("hrbrthemes")
options(tibble.width = Inf)
library("Matrix")
#install.packages("Seurat")
library("Seurat")


set.seed(2)
t1 <- Sys.time()
rds_dir <- "rds"
dir.create(rds_dir, showWarnings = FALSE)


# load data ===============

message("Reading metadata...")
meta <- read.csv("2093-Olbrecht_metadata.csv")

message("Loading matrix...")
z <- readRDS(file = "2095-Olbrecht_counts_matrix.rds")

message("Metadata per cell...")
cell_names <- data.frame(cell_label = colnames(z),
                         sample_name = sub(".*_(.*)",
                                           "\\1",
                                           colnames(z)),
                         stringsAsFactors = FALSE,
                         row.names = colnames(z))
stopifnot(all(cell_names$sample_name %in% meta$sample_name))
meta_seurat <- left_join(x = cell_names,
                         y = as.data.frame(meta),
                         by = "sample_name")
# Seurat needs rownames, add them manually
# left_join preserves the order of cells:
stopifnot(all.equal(meta_seurat$cell_label, cell_names$cell_label))
# add cell labels as rownames
rownames(meta_seurat) <- meta_seurat$cell_label

# Seurat pipeline ===============

message("Creating Seurat object...")
x <- CreateSeuratObject(count = z,
                        project              = "OV",
                        min.cells            = 10,     # only genes > 10 cells
                        min.genes            = 200,    # only cells with > 200 genes
                        is.expr              = 0,      # expression threshold for 'detected' gene
                        normalization.method = NULL,   # no normalization for now
                        scale.factor         = 10000,  # scale factor in normalization process
                        # (not used given norm.method == NULL)
                        do.scale             = FALSE,  # gene-based z-score
                        do.center            = FALSE,  # gene-based centering
                        names.field          = 1,
                        names.delim          = "_",
                        meta.data            = meta_seurat,
                        display.progress     = TRUE)

# save data ===============

message()
message("Saving data...")
saveRDS(x, file = file.path(rds_dir, "rds_12_data_seurat.rds"))






message()
message("Done.\n")
print(Sys.time() - t1)
