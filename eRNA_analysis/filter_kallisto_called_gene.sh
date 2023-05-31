#!/bin/bash

#sort -k1,1 -k2,2n ph2_enh.bed > ph2_enh.sort.bed

#bedtools intersect -wa -wb \
# -a ph2_enh.sort.bed \
# -b gene_hg38.sort.bed \
# -sorted \
# -v > ph2_enh_gene.sort.bed

#sort -k1,1 -k2,2n ph1_random.bed > ph1_random.sort.bed
#sort -k1,1 -k2,2n ph2_random.bed > ph2_random.sort.bed

bedtools intersect -wa -wb \
 -a ph1_random.sort.bed \
 -b gene_hg38.sort.protein.coding.bed \
 -sorted \
 -v > ph1_random_prot_coding.sort.bed

bedtools intersect -wa -wb \
 -a ph2_random.sort.bed \
 -b gene_hg38.sort.protein.coding.bed \
 -sorted \
 -v > ph2_random_prot_coding.sort.bed

bedtools intersect -wa -wb \
-a ph1_enh.sort.bed \
-b gene_hg38.sort.protein.coding.bed \
-sorted \
-v > ph1_enh_prot_coding.sort.bed

bedtools intersect -wa -wb \
-a ph2_enh.sort.bed \
-b gene_hg38.sort.protein.coding.bed \
-sorted \
-v > ph2_enh_prot_coding.sort.bed
