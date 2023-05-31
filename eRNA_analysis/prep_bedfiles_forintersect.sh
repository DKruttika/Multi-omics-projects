# get the enhancer intervals in both POCROC and pancan data sets

# for the pancan data, took the prognostic enhancers found in that paper

# used bedtools intersect by installing bedtools in miniconda

# tab separated .bed file
awk -v OFS="\t" '$1=$1' pancan_prognostic_eRNA.bed > pancan_prognostic_eRNA_tab.bed
awk -v OFS="\t" '$1=$1' sign_enh_calc_25k_64k_50923.bed > sign_enh_calc_25k_64k_50923_tab.bed

# sorted bed file; sorting needs to be chr1, chr10, chr11, and not 1,2,3

awk 'NR>1{print $0}' pancan_prognostic_eRNA_tab.bed | sort -k1,1 -k2,2n > pancan_prognostic_eRNA_tab_sort.bed
awk 'NR>1{print $0}' sign_enh_calc_25k_64k_50923_tab.bed | sort -k1,1 -k2,2n > sign_enh_calc_25k_64k_50923_tab_sort.bed
