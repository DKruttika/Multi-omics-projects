############################################################
##              reading the data set                     ##
############################################################
rm(list = ls(all.names = TRUE))
gc()

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))


cell_type_col <-  c("B cells" ="#0000FF", "Cytotoxic T cells" ="#FF0000", "Endothelial cells" ="#00FF00",
                    "Fibroblasts" ="#000033","Helper T cells"= "#FF00B6", "M1 Macrophages"= "#005300",
                    "M2 Macrophages" ="#FFD300","Necrotic Tumor cells"= "#009FFF", "other"="#9A4D42",
                    "Smooth muscle cells"="#00FFBE", "Tregs"="#783FC1","Tumor cells"= "#1F9698",
                    "Dendritic cells"="#FFACFD", "NK cells"="#B1CC71","Monocytes"= "#F1085C",
                    "T cells double negative"="#FE8F42", "Epithelial cells"="#DD00FF")

base_dir <- "/path/box/folder"

args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]

start_time <- Sys.time()
df <- read.csv(paste0(base_dir, "AKOYA_Data_ROI_specific_92223/", fn)) 

filename <- fn

str1 <- strsplit(filename, "_")[[1]][7]
str2 <- strsplit(filename, "_")[[1]][8]
str3 <- strsplit(filename, "_")[[1]][9]
paste1 <- paste(str1, str2, sep = "-")
paste2 <- paste(paste1, str3, sep = "_")

metadata <- read_csv(paste0(base_dir, "Akoya_sample_info_25slides_072023.csv"))
metadata <- as.data.frame(metadata)
# S13-32894-G1 is labelled as S13-32894_G1 in the merged file
# changing this specific accession number in the metadata file

metadata[grep("S13-32894-G1", metadata$accession_number),1] <- "S13-32894_G1"

output_filename = (filter(metadata, accession_number == paste2))$POCROC_labels

png(paste0("./Figures/tSNE_plots/", output_filename, ".png"), units = "in", height = 8, width = 12, res = 300)
ggplot(df, aes(x=tSNE1, y=tSNE2, col = df[,ncol(df)])) +
  geom_point(size=2) + 
  theme_minimal() +
  scale_color_manual(values = cell_type_col) +
  labs(title=paste((filter(metadata, accession_number == paste2))$image_svs_name,
                   (filter(metadata, accession_number == paste2))$POCROC_labels, sep = "/"), 
       x="t-SNE 1", y="t-SNE 2")
dev.off()  

end_time <- Sys.time()
time_taken <- round(end_time - start_time,2)
print(time_taken)


