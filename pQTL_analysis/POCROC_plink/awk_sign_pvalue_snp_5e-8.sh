#!/bin/bash



awk -F'\t' -v OFS='\t' -v COL_NUM=$1 -v PROT_NAME=$2 '{if($COL_NUM < 0.00000005) {print $1, $2, $3, $4, PROT_NAME}}' ./pvalues_pqtl/combined_sign_plink_pvalue_5e-8_92021.txt >> temp_pocroc_sign_snp_5e-8.txt
