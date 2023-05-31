#!/bin/bash


#grep -w "PAX8" GeneHancer_AnnotSV_gene_association_scores_v5.15.txt | awk '{print $2}' | sort | uniq -c
#awk '{print $4}' | sort | uniq -c

grep -w "MYC" GeneHancer_AnnotSV_gene_association_scores_v5.15.txt | grep -v "lnc" | awk '{if ($4==1) {print $0}}' > MYC_genehancer.txt
grep -w "PAX8" GeneHancer_AnnotSV_gene_association_scores_v5.15.txt | grep -v -E "AS1|lnc" | awk '{if ($4==1) {print $0}}' > PAX8_genehancer.txt
grep -w "KRAS" GeneHancer_AnnotSV_gene_association_scores_v5.15.txt | grep -v "lnc" | awk '{if ($4==1) {print $0}}' > KRAS_genehancer.txt
grep -w "MUC16" GeneHancer_AnnotSV_gene_association_scores_v5.15.txt | grep -v "lnc" | awk '{if ($4==1) {print $0}}' > MUC16_genehancer.txt
grep -w "WT1" GeneHancer_AnnotSV_gene_association_scores_v5.15.txt | grep -v "AS" | awk '{if ($4==1) {print $0}}' > WT1_genehancer.txt

