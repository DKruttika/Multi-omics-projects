# AKOYA_Data_ROI_specific_92223
rm(list = ls(all.names = TRUE))
gc()


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tibble))

base_dir <- "/path/box/folder"

args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]

start_time <- Sys.time()
df <- read.csv(paste0(base_dir, fn))
end_time <- Sys.time()

time_taken <- round(end_time - start_time,2)
print(time_taken)

print(table(df$ROI_ID))

df <- df[grep("Tumor", df$ROI_ID),] 

df_table <- as.data.frame(table(df[,ncol(df)]))
df_table$count_norm <- df_table$Freq/nrow(df)
filename <- fn

str1 <- strsplit(filename, "_")[[1]][7]
str2 <- strsplit(filename, "_")[[1]][8]
str3 <- strsplit(filename, "_")[[1]][9]

paste1 <- paste(str1, str2, sep = "-")
paste2 <- paste(paste1, str3, sep = "_")

df_table$accession_number <- paste2

results_dir <- "/path/box/folder"


write.table(df_table, paste0(results_dir, fn),
            sep = "\t", row.names = F, col.names = T, quote = F)

print(paste0("Generated summarized results for ", paste2))

