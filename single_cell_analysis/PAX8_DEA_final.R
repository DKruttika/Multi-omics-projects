#Fresh biopsies were obtained from these 7 treatment-naive patients during primary cytoreductive 
#surgery (patients P1 and P4) or during diagnostic laparoscopy (patients P2–P3, P5–P7) and consisted 
#of primary ovarian tumour (P1, P4), intraperitoneal metastatic lesions (peritoneum (P1–P3, P5–P7) or 
#omentum (P1)) and normal adjacent tissue (P1 (omental and peritoneal tissue), P4 (ovary))

#except patient P3 who presented with a mixed ovarian epithelial carcinoma consisting of 
#clear cell and high-grade serous components


##### identify tumor cell clusters

x.markers <- FindAllMarkers(x, only.pos = TRUE, min.pct = 0.25)
top10_markers <- x.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

table(x@active.ident)

library(readxl)
regev_malignant <- read_excel("regev_malignant_markers.xlsx")

grep("KLK", regev_malignant$`Cell type`)

(regev_malignant[grep("KLK", regev_malignant$`Cell type`),])$`Cell type`
(regev_malignant[grep("KRT", regev_malignant$`Cell type`),])$`Cell type`

tumor_clusters <- x.markers %>% dplyr::group_by(cluster) %>% dplyr::filter(., gene %in% c("EPCAM","CD24","PAX8",(regev_malignant[grep("KLK", regev_malignant$`Cell type`),])$`Cell type`,
                                                                                (regev_malignant[grep("KRT", regev_malignant$`Cell type`),])$`Cell type`)) %>% 
  filter(., p_val_adj < 0.05) %>% filter(., avg_log2FC > 0.2)

table(x$sample_type, x$seurat_clusters)

unique((filter(x.markers, p_val_adj < 0.05) %>% filter(., gene %in% (regev_malignant[grep("KLK", regev_malignant$`Cell type`),])$`Cell type` ))$gene)


#3, 6, 8, 12, 14, 19, 22, 23 - clusters that had mainly EPCAM expression along with KRT genes
# and all clusters except 3 have KLK genes overexpressed along withe EPCAM and KRT
# given first priority to EPCAM and then KRT genes followed by KLK

write.csv(x.markers, "cluster_DEA.csv")

# plotting different gene expressions on tsne

cluster <- readRDS("rds/rds_15_data_clusters.rds")



# expressed everywhere: KRT10, KRTCAP2
FeaturePlot(x, features = c("KRT18"))

new.cluster.ids <- c("random","random", "random", "Tumor3", "random","Tumor5","Tumor6", "random",
                     "Tumor8","random","random","random","Tumor12", "random", "Tumor14",
                     "random","random","random","random","Tumor19", "random","random","Tumor22",
                     "Tumor23")
names(new.cluster.ids) <- levels(cluster)
cluster <- RenameIdents(cluster, new.cluster.ids)


kelly.colours <- c("gray95", "gray13", "gold2", "plum4", "darkorange1", "lightskyblue2", "firebrick", "burlywood3", "gray51", "springgreen4", "lightpink2", "deepskyblue4", "lightsalmon2", 
                   "mediumpurple4", "orange", "maroon", "yellow3", "brown4", "yellow4", "sienna4", "chocolate", "gray19")

png("TSNE_tumor_cells.png",units = "in", height = 6, width = 6, res = 600)
DimPlot(object = cluster, cols = kelly.colours[1:10], pt.size = 0.5)
dev.off()

png("TSNE_features.png",units = "in", height = 10, width = 12, res = 600)
FeaturePlot(cluster, features = c("EPCAM", "KRT18", "CD24", "KLK6", "PAX8"))
dev.off()

cluster <- readRDS("rds/rds_15_data_clusters.rds")
cluster <- SetIdent(cluster, value = cluster$patient_id)

png("TSNE_patientID_cells.png",units = "in", height = 6, width = 6, res = 600)
DimPlot(object = cluster, cols = c('#9e0142','#000000','#66c2a5','#5e4fa2','#fb9a99','#a6cee3','#ff7f00'), pt.size = 0.5)
dev.off()

########## take only the define tumor clusters
# also remove the 9 normal cells
table(cluster$sample_type, cluster$seurat_clusters)


DEA_df <- subset(cluster, sample_type == "Tumor")
table(DEA_df$sample_type)

DEA_df <- subset(DEA_df, subset = seurat_clusters %in% c(3,5, 6, 8, 12, 14, 19, 22, 23)) # these clusters have 11 cells from "normal" samples
table(DEA_df$sample_type)

DEA_df <- SetIdent(DEA_df, value = ifelse(DEA_df@assays$RNA@data["PAX8", ] > 0, "Present", "Absent"))
DEA_df$pax8_levels <- ifelse(DEA_df@assays$RNA@data["PAX8", ] > 0, "Present", "Absent")

table(DEA_df$pax8_levels) # 2848 cells dont have PAX8 expression (46%) out of 6116
table(DEA_df$patient_id, DEA_df$pax8_levels)
table(DEA_df$seurat_clusters, DEA_df$pax8_levels)

PAX8.markers <- FindMarkers(DEA_df, ident.1 = "Absent", ident.2 = "Present")
write.csv(PAX8.markers, "PAX8_DEA.csv")



write.table(rownames(filter(PAX8.markers, p_val_adj < 0.05) %>% filter(., avg_log2FC < -0.5)),
            "PAX8_negLFC_regulation_single_cell.txt", sep = "\t",quote = F, col.names = F,
            row.names = F)

write.table(rownames(filter(PAX8.markers, p_val_adj < 0.05) %>% filter(., avg_log2FC > 0.5)),
            "PAX8_posLFC_single_cell.txt", sep = "\t",quote = F, col.names = F,
            row.names = F)



