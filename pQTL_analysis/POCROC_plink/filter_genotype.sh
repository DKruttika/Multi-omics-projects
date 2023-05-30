#!/bin/bash

# input should be:
# 1: column number of protein
# 2: protein name
# second script
# 3: chromosome_snp
# 4: protein name 

OUTPUT="$2_genotype.txt"

awk -F'\t' -v COL_NUMBER=$1 -v PROT_NAME=$2 '{if(($COL_NUMBER <= 0.00000001)){print "chr" $1"_"$3" "PROT_NAME}}' ./pvalues_pqtl/combined_sign_plink_pvalue2_filt_91721.txt >> ${OUTPUT}

cat *_genotype.txt | awk -v OFS='\t' '{ print $1, $2 }' > sign_snp_location_plink_POCROC.txt
 
# will have to run this separately
#OUTPUT2="POCROC_sign_genotype.txt"

#awk -F'\t' -v chr_num=$1 -v PROT_NAME=$2 '{if($1 == chr_num) {print $0"\t"PROT_NAME}}' genotype_26samples_ph1_ph2_60921.txt >> ${OUTPUT2} 

