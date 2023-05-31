rm(list=ls())

library("tidyverse")
#library("scales")
#library("gridExtra")
#library("hrbrthemes")
#options(tibble.width = Inf)
library("Seurat")

t1 <- Sys.time()
set.seed(2)
rds_dir <- "rds"

# load data ===============

message("Loading data...")
x <- readRDS(file = file.path(rds_dir, "rds_12_data_seurat.rds"))

# mitochondrial genes ===============

message("Evaluating mitochondrial genes...")
mito.genes <- grep("^MT-", rownames(x@assays$RNA), value = TRUE)
print(paste(sort(mito.genes), collapse = ", "))
message("found ", length(mito.genes), " mitochondrial genes\n", sep = "")

percent.mito <- Matrix::colSums(x@assays$RNA[mito.genes]) /
  Matrix::colSums(x@assays$RNA)
stopifnot(all.equal(names(percent.mito), colnames(x@assays$RNA)))

x <- AddMetaData(object = x,
                 metadata = percent.mito,
                 col.name = "percent.mito")
# filter cells ===============

# Filter out cells that have unique gene counts over 6,000 or mitochondrial
# content over 15%
message("Subsetting data...")

x <- subset(x = x, subset = nFeature_RNA < 6000 & percent.mito < 0.15)

# normalize data ===============

x <- NormalizeData(object = x,
                   normalization.method = "LogNormalize",
                   scale.factor = 10000)
# variable genes ===============

message("Variable genes...")
x_low_cutoff  <- 0.0125
x_high_cutoff <- 3 # (~ outliers)
y_cutoff      <- 0.5  # ----------------> z-score > 0.5
y_high_cutoff <- Inf
x <- FindVariableFeatures(object = x,
                          mean.function = ExpMean,
                          dispersion.function = LogVMR,
                          mean.cutoff = c(x_low_cutoff, x_high_cutoff),
                          dispersion.cutoff = c(y_cutoff, y_high_cutoff)) 

# cell cycle ===============

message("Scores for cell cycle...")
# http://satijalab.org/seurat/cell_cycle_vignette.html:
# Read in a list of cell cycle markers, from Tirosh et al, 2015.
# See Seurat website:
cc.genes <- readLines(con = "regev_lab_cell_cycle_genes.txt")

# ssegregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
s.genes <- ifelse(s.genes == "MLF1IP", "CENPU", s.genes)

g2m.genes <- cc.genes[44:97]


x <- CellCycleScoring(object = x, s.features = s.genes,
                      g2m.features = g2m.genes,
                      set.ident = TRUE)
#head(y = y@meta.data)

# scaling data and removing sources of unwanted variation ===============

message("Scaling data and removing sources of unwanted variation...")
x <- ScaleData(x, vars.to.regress = c("nCount_RNA",
                                      "percent.mito",
                                      "S.Score",
                                      "G2M.Score"))

# save data ===============

message("Saving data...")
saveRDS(x, file = file.path(rds_dir, "rds_13_data_filtered.rds"))






message("Done.")
print(Sys.time() - t1)





