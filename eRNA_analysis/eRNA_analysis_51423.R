ph1_col_header <- read.table("./eRNA_data/ph1_column_header.txt", sep = "\t")
ph1_col_header <- separate(ph1_col_header, col = V1, into = "sample", sep = "-")

ph2_col_header <- read.table("./eRNA_data/ph2_column_header.txt", sep = "\t")
ph2_col_header <- separate(ph2_col_header, col = V1, into = c("sample", "sample2"), sep = "-")
ph2_col_header$sample <- paste(ph2_col_header$sample, ph2_col_header$sample2, sep = ".")


ph1_enh <- read.table("./eRNA_data/ph1_enh_prot_coding.sort.bed", sep = "\t", header = F)
ph1_enh$enhancer <- paste(ph1_enh$V1, ph1_enh$V2, sep = ":")
ph1_enh$enhancer <- paste(ph1_enh$enhancer, ph1_enh$V3, sep = "-")
ph1_enh <- subset(ph1_enh, select = -c(V1, V2, V3))
colnames(ph1_enh) <- c(ph1_col_header$sample, "enhancer")
ph1_enh <- column_to_rownames(ph1_enh, "enhancer")

ph2_enh <- read.table("./eRNA_data/ph2_enh_prot_coding.sort.bed", sep = "\t", header = F)
ph2_enh$enhancer <- paste(ph2_enh$V1, ph2_enh$V2, sep = ":")
ph2_enh$enhancer <- paste(ph2_enh$enhancer, ph2_enh$V3, sep = "-")
ph2_enh <- subset(ph2_enh, select = -c(V1, V2, V3))
colnames(ph2_enh) <- c(ph2_col_header$sample, "enhancer")
ph2_enh <- column_to_rownames(ph2_enh, "enhancer")

all(ph1_enh$enhancer == ph2_enh$enhancer)
ph1_ph2_enh <- cbind(ph1_enh, ph2_enh)
ph1_ph2_enh <- ph1_ph2_enh[,good.samples_ph1_ph2$sample]

ph1_ph2_enh_filt <- ph1_ph2_enh[apply(ph1_ph2_enh == 0, 1, sum) <= 49, ]




ph1_random <- read.table("./eRNA_data/ph1_random_prot_coding.sort.bed", sep = "\t", header = F)
ph1_random$enhancer <- paste(ph1_random$V1, ph1_random$V2, sep = ":")
ph1_random$enhancer <- paste(ph1_random$enhancer, ph1_random$V3, sep = "-")
ph1_random <- subset(ph1_random, select = -c(V1, V2, V3))
colnames(ph1_random) <- c(ph1_col_header$sample, "enhancer")
ph1_random <- column_to_rownames(ph1_random, "enhancer")

ph2_random <- read.table("./eRNA_data/ph2_random_prot_coding.sort.bed", sep = "\t", header = F)
ph2_random$enhancer <- paste(ph2_random$V1, ph2_random$V2, sep = ":")
ph2_random$enhancer <- paste(ph2_random$enhancer, ph2_random$V3, sep = "-")
ph2_random <- subset(ph2_random, select = -c(V1, V2, V3))
colnames(ph2_random) <- c(ph2_col_header$sample, "enhancer")
ph2_random <- column_to_rownames(ph2_random, "enhancer")

all(rownames(ph1_random) == rownames(ph2_random))

ph1_ph2_random <- cbind(ph1_random, ph2_random)
ph1_ph2_random <- ph1_ph2_random[,good.samples_ph1_ph2$sample]

ph1_ph2_random_filt <- ph1_ph2_random[apply(ph1_ph2_random == 0, 1, sum) <= 49, ]


all(colnames(ph1_ph2_enh) == colnames(ph1_ph2_random))

# calculating significant enhancers:

library(edgeR)
scale_ph1_ph2_enh <- calcNormFactors(ph1_ph2_enh_filt) 
scale_ph1_ph2_enh <- sweep(as.matrix(ph1_ph2_enh_filt), 2, scale_ph1_ph2_enh, "/")
scale_ph1_ph2_enh <- as.data.frame(scale_ph1_ph2_enh)

scale_random <- calcNormFactors(ph1_ph2_random_filt) 
scale_random <- sweep(as.matrix(ph1_ph2_random_filt), 2, scale_random, "/")
scale_random <- as.data.frame(scale_random)

df_enh <- tibble()

for (j in 1:50) {
  for (i in 1:nrow(scale_ph1_ph2_enh)) {
    df_enh[i,j] <- sum(scale_random[,j] > scale_ph1_ph2_enh[,j][i])
    #significance <- result_enh_2[i]
  }
}

colnames(df_enh) <- colnames(scale_ph1_ph2_enh)
rownames(df_enh) <- rownames(scale_ph1_ph2_enh)

df_enh <- df_enh/nrow(scale_random)

sign_pvalue <- apply(df_enh < 0.01, 1, sum)
sign_pvalue <- as.data.frame(sign_pvalue)
rownames(sign_pvalue) <- NULL
sign_pvalue$enh <- rownames(scale_ph1_ph2_enh)
colnames(sign_pvalue) <- c("sign", "enh")

write.csv(df_enh, "./eRNA_data/sign_enh_calc_25k_64k_50923.csv")

# only 275 enhancers significant across all samples at pvalue < 0.01

sign_pvalue_50 <- filter(sign_pvalue, sign == 50)



############################################################
##      barplot of significant enhancers for every sample            
############################################################
sign_enh_count <- as.data.frame(colSums(df_enh < 0.01))
colnames(sign_enh_count) <- "enh_count"

sign_enh_count <- rownames_to_column(sign_enh_count, "sample")
sign_enh_count <- left_join(sign_enh_count, good.samples_ph1_ph2)

sign_enh_count <- sign_enh_count[order(sign_enh_count$patientID),]
sign_enh_count <- mutate(sign_enh_count, patientID = as.factor(patientID))

sign_enh_count$prim_HRD <- paste(sign_enh_count$isRecur, sign_enh_count$HRD_status, sep = "-")

sign_enh_count$HRD_num <- ifelse(sign_enh_count$HRD_status == "HRD", 0, 1)

png("./figures/count_sign_enh_per_sample.png", units = "in", height = 12, 
    width = 22, res = 600)
ggplot(sign_enh_count, aes(x = forcats::fct_reorder(patientID, HRD_num), y = enh_count, 
                           fill = as.character(recurNumber))) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#f97c03","#9E95C7","#7E72B4","#5e4fa2", "#2E2751")) +
  theme_bw() +
  theme(text = element_text(size = 25)) +
  labs(y = "count(Significant eRNAs)", x = "Patient IDs", fill = "Recurrence")
#scale_fill_manual(values = c("#AA323F",'#286D97', "#E99DA5",'#96C2DD'))
dev.off()

############################################################
##      boxplots of counts of eRNAs per 4 groups
############################################################

my_comparisons <- list( c("Primary-HRD", "Recurrance-HRD"),
                        c("Primary-HRP", "Recurrance-HRP"),
                        c("Primary-HRD", "Primary-HRP"),
                        c("Recurrance-HRD", "Recurrance-HRP"))
comparisons_HRD <- list( c("HRD", "HRP"))
comparisons_prim <- list( c("Primary", "Recurrance"))

sign_enh_count <- sign_enh_count[order(sign_enh_count$prim_HRD),]
sign_enh_count <- mutate(sign_enh_count, prim_HRD = as.factor(prim_HRD))


library(ggpubr)
p <- ggboxplot(sign_enh_count, x = "prim_HRD", y = "enh_count",
               fill = "prim_HRD", palette = c("#AA323F",'#286D97', "#E99DA5",'#96C2DD'),
               add = "jitter",
               ylab = "count(Significant eRNAs)")

####### for recurrence only
sign_enh_count <- sign_enh_count[order(sign_enh_count$isRecur),]
sign_enh_count <- mutate(sign_enh_count, isRecur = as.factor(isRecur))


png("./figures/count_sign_enh_recur.png", units = "in", height = 8, 
    width = 8, res = 600)
p + stat_compare_means(comparisons = c(my_comparisons), 
                       method = "t.test",  size = 8) +
  theme(text = element_text(size = 25),
        axis.text.x=element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top") +
  labs(fill = "")
dev.off()

############################################################
##      boxplots of expression of eRNAs per 4 groups individual enhancers
############################################################
ph1_ph2_enh_sign <- filter(ph1_ph2_enh_filt, rownames(ph1_ph2_enh_filt) %in% 
                               sign_pvalue_50$enh)
ph1_ph2_enh_sign <- rownames_to_column(ph1_ph2_enh_sign, "eRNAs")

ph1_ph2_enh_sign <- gather(ph1_ph2_enh_sign, key = "sample", value = "eRNA_exp", -eRNAs)
ph1_ph2_enh_sign <- left_join(ph1_ph2_enh_sign, good.samples_ph1_ph2)

ph1_ph2_enh_sign$prim_HRD <- paste(ph1_ph2_enh_sign$isRecur, ph1_ph2_enh_sign$HRD_status,
                                     sep = "-")

ph1_ph2_enh_sign <- ph1_ph2_enh_sign[order(ph1_ph2_enh_sign$prim_HRD),]
ph1_ph2_enh_sign <- mutate(ph1_ph2_enh_sign, prim_HRD = as.factor(prim_HRD))

p <- ggboxplot(ph1_ph2_enh_sign, x = "prim_HRD", y = "eRNA_exp",
               fill = "prim_HRD", palette = c("#AA323F",'#286D97', "#E99DA5",'#96C2DD'),
               add = "jitter",
               ylab = "TPM(Significant eRNAs)")

png("./figures/exp_each_sign_enh.png", units = "in", height = 10, 
    width = 12, res = 600)
p + stat_compare_means(comparisons = c(my_comparisons), 
                       method = "t.test",  size = 8) +
  theme(text = element_text(size = 25),
        axis.text.x=element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top") +
  labs(fill = "")
dev.off()

############################################################
##      heatmap of significant eRNAs across all samples            
############################################################
ph1_ph2_enh_sign <- filter(ph1_ph2_enh_filt, rownames(ph1_ph2_enh_filt) %in% 
                             sign_pvalue_50$enh)
all(colnames(ph1_ph2_enh_sign) == rownames(metadata_rna))

# metadata and complexheatmap from line 45 of this script

ph1_ph2_enh_sign <- as.data.frame(t(scale(t(ph1_ph2_enh_sign))))

hmap <- Heatmap(
  as.matrix(ph1_ph2_enh_sign),
  #name = "Protein Intensity(log2)",
  show_row_names = F,
  show_column_names = F,
  cluster_rows = T,
  cluster_columns = T,
  show_column_dend = TRUE,
  show_row_dend = T,
  row_dend_reorder = T,
  column_dend_reorder = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  width = unit(100, "mm"), #col = Grays,
  top_annotation=colAnn,
  heatmap_legend_param = list(
    title = "scaled(sign eRNA expression)",
    title_position = "lefttop-rot"))

png("eRNA_all_sign_pval_0.01.png",height=8,width=12, res = 600,units = "in") # change height and width parameter
ht <- draw(hmap, heatmap_legend_side="right", annotation_legend_side="left")
dev.off()


# number of enhancers detected in more than 80% of samples
# 6060 enhancers quantified across 80% samples

dim(ph1_ph2_enh_filt[apply(ph1_ph2_enh_filt == 0, 1, sum) <= 10, ])

############################################################
# prognostic enhancers identified in the pancan paper
############################################################

pancan_common <- read.table("./eRNA_data/filt_common_pancan_POCROC_sort.bed", sep = "\t")

pancan_common$pancan_enhancer <- paste(pancan_common$V1, pancan_common$V2, sep = ":")
pancan_common$pancan_enhancer <- paste(pancan_common$pancan_enhancer, pancan_common$V3, sep = "-")

pancan_common$pocroc_enhancer <- paste(pancan_common$V4, pancan_common$V5, sep = ":")
pancan_common$pocroc_enhancer <- paste(pancan_common$pocroc_enhancer, pancan_common$V6, sep = "-")



pancan_prognostic <- read.table("./eRNA_data/pancan_prognostic_eRNA.txt", sep = "\t", header = T)


intersect(pancan_common$pancan_enhancer, pancan_prognostic$Prognostic_enhancer)

sign_enh <- read.csv("./eRNA_data/sign_enh_calc_25k_64k_50923.csv")


sign_pvalue <- apply(sign_enh < 0.01, 1, sum)
sign_pvalue <- as.data.frame(sign_pvalue)
rownames(sign_pvalue) <- NULL
sign_pvalue$enh <- rownames(ph1_ph2_enh_filt)
colnames(sign_pvalue) <- c("sign", "enh")

intersect(sign_pvalue_50$enh, pancan_common$pocroc_enhancer) #1 chr11:341474-341670

intersect(unique(pancan_common$pancan_enhancer), pancan_prognostic$Prognostic_enhancer) #269

table((filter(pancan_prognostic, Prognostic_enhancer %in% unique(pancan_common$pancan_enhancer)))$Disease)

# there is an overlap between ovarian cancer prognostic eRNA and the one identified in POCROC data
# overlap is of 278 bp out of the 390 in each enhancer region; 71%

pancan_disease_overlap <- as.data.frame(table((filter(pancan_prognostic, Prognostic_enhancer %in% unique(pancan_common$pancan_enhancer)))$Disease))

pancan_disease_overlap <- pancan_disease_overlap[order(pancan_disease_overlap$Freq, decreasing = T),]

png("common_pancan_prognostic.png",height=8,width=12, res = 600,units = "in") # change height and width parameter
ggplot(pancan_disease_overlap, aes(y = Freq, 
                                   x = forcats::fct_reorder(Var1, Freq, .desc = TRUE), fill = Var1)) +
  geom_col() +
  theme_bw() +
  labs(x = "Prognostic enhancers", y = "Overlapping enhancers") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none",
        text = element_text(size = 15)) +
  scale_fill_manual(values = c('#FAF7F0', '#BDD5DD', '#2C2D73', '#737169', 
                               '#2C3823', '#769F3C', '#D8B3C2', '#E8E451', 
                               '#A8A459', '#A1BC73', '#E9D2B1', '#532F73', 
                               '#B38EB7', '#B6AB83', '#64221E', '#861F61', 
                               '#21215B', '#E1EADA', '#436669', '#E4CE48', 
                               '#4F1715', '#49397F'))

dev.off()

############################################################
#### known enhancers - Genehancer data
############################################################


specific_enh <- read.table("./eRNA_data/filt_common_genehancer_enh_POCROC_sort.bed", sep = "\t")
specific_enh$specific_enhancer <- paste(specific_enh$V1, specific_enh$V2, sep = ":")
specific_enh$specific_enhancer <- paste(specific_enh$specific_enhancer, specific_enh$V3, sep = "-")

specific_enh$pocroc_enhancer <- paste(specific_enh$V4, specific_enh$V5, sep = ":")
specific_enh$pocroc_enhancer <- paste(specific_enh$pocroc_enhancer, specific_enh$V6, sep = "-")

length(unique(specific_enh$pocroc_enhancer))

intersect(sign_pvalue_50$enh, specific_enh$pocroc_enhancer) 

specific_eRNA <- filter(ph1_ph2_enh_filt, rownames(ph1_ph2_enh_filt) %in% 
                               specific_enh$pocroc_enhancer)

all(colnames(specific_eRNA) == rownames(metadata_rna))

gene_enh <- as.data.frame(specific_enh$pocroc_enhancer)
colnames(gene_enh) <- "enhancers"
gene_enh$gene <- "PAX8"
gene_enh$gene[grep("chr8", gene_enh$enhancers)] <- "MYC"
gene_enh$gene[grep("chr19", gene_enh$enhancers)] <- "MUC16"
#gene_enh$gene[grep("chr12", gene_enh$enhancers)] <- "KRAS"

gene_enh <- column_to_rownames(gene_enh, "enhancers")

all(rownames(gene_enh) == rownames(specific_eRNA))
condition_colors2 <- list(gene = c("#5e4fa2", "#66c2a5", "#fdae61"))
names(condition_colors2$gene) <- unique(gene_enh$gene)

colAnn2 <- HeatmapAnnotation(df=gene_enh, which="row",
                             col = condition_colors2,
                             show_annotation_name = FALSE,
                             annotation_width=unit(c(2, 8), "cm"), 
                             gap=unit(3, "mm"), annotation_name_gp= gpar(fontsize = 15))


specific_eRNA_scale <- as.data.frame(t(scale(t(specific_eRNA))))

lajolla <- sequential_hcl(1000, "Lajolla")

library(viridis)

hmap <- Heatmap(
  as.matrix(specific_eRNA),
  #name = "Protein Intensity(log2)",
  show_row_names = F,
  show_column_names = F,
  cluster_rows = T,
  cluster_columns = T,
  show_column_dend = TRUE,
  show_row_dend = T,
  row_dend_reorder = T,
  column_dend_reorder = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  width = unit(100, "mm"), col = viridis::cividis(500),
  top_annotation=colAnn, right_annotation = colAnn2,
  heatmap_legend_param = list(
    title = "tpm(eRNA expression)",
    title_position = "lefttop-rot"))

png("eRNA_MYC_PAX8.png",height=8,width=12, res = 600,units = "in") # change height and width parameter
ht <- draw(hmap, heatmap_legend_side="right", annotation_legend_side="left")
dev.off()


sign_enh_275 <- as.data.frame(sign_pvalue_50$enh)
colnames(sign_enh_275) <- "enh"

sign_enh_275 <- separate(sign_enh_275, col = enh, into = c("chr", "random"), sep = ":")
sign_enh_275 <- separate(sign_enh_275, col = random, into = c("start", "end"), sep = "-")

write.table(sign_enh_275, "sign_enh_275.bed", sep = "\t", col.names = F, row.names = F, quote = F)



