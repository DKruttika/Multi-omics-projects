#!/bin/bash

awk '{if($3=="exon"){print $0}}' gencode.v34.annotation.gtf > exon_hg38.gtf
#awk '{if($3=="transcript"){print $0}}' gencode.v34.annotation.gtf > transcript_hg38.gtf

awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' exon_hg38.gtf | gtf2bed -> exon_hg38.bed

#awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' transcript_hg38.gtf | gtf2bed -> t

# downloaded the entire f5 enhancer list as is file name: enhancer_fantom5_hg38_og.bed
# got this file from R after separating the columns into chrm no. start and end

sort -k1,1 -k2,2n F5.hg38.enhancers.bed > F5.hg38.enhancers.sort.bed

#sort -k1,1 -k2,2n hg38-blacklist.v2.bed > hg38-blacklist.v2_sort.bed

sort -k1,1 -k2,2n bed_file_kallisto_hg38.bed > pancan_enh_hg38_sort.bed

bedtools intersect -wa -wb \
 -a F5.hg38.enhancers.sort.bed \
 -b exon_hg38.bed hg38-blacklist.v2_sort.bed \
 -sorted \
 -v > exon.gencode.F5.hg38.enhancers.sort.bed

bedtools intersect -wa -wb \
 -a pancan_enh_hg38_sort.bed \
 -b exon_hg38.bed hg38-blacklist.v2_sort.bed \
 -sorted \
 -v > exon.gencode.pancan.bed


# to get the unique regions in the pancan list which are not in my filtered list of enhancers
bedtools intersect -v -a pancan_enh_hg38_sort.bed -b exon.gencode.F5.hg38.enhancers.sort.bed > test2.bed

# merging the two files to get a complete list of enhancers

cat exon.gencode.F5.hg38.enhancers.sort.bed test2.bed > merge.F5.hg38.pancan.bed

sort -k1,1 -k2,2n merge.F5.hg38.pancan.bed > merge.F5.hg38.pancan.sort.bed

awk '{if($3-$2 >= 151) print}' merge.F5.hg38.pancan.sort.bed > merge.F5.hg38.pancan.sort.filt150.bed  

