#!/bin/bash


awk -F'\t' -v OFS='\t' -v COL_NUM=$1 -v PROT_NAME=$2 '{if($COL_NUM < 0.00000005) {print $1, $2, $3, $COL_NUM, PROT_NAME}}' combined_sign_TCGA.txt >> temp_tcga_sign_snp_5e-8_pvalue.txt
