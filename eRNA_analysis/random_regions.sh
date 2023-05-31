#!/bin/bash

# use bedtools shuffle to create random regions of the genome
# use triple that of the 53K regions as input

# filter those regions for exon, blacklist and eRNA regions as well as intervals of repeats in the genome

# input for shuffles: merge.F5.hg38.pancan.sort.filt150.interval.bed times 3

# input for exon list: exon_hg38.bed

# input for blacklist regions: hg38-blacklist.v2_sort.bed

# input for repeat intervals: 

#cat merge.F5.hg38.pancan.sort.filt150.bed | awk 'BEGIN{OFS="\t";} {print $1,$2,$3;}' > merge.F5.hg38.pancan.sort.filt150.interval.bed

cat merge.F5.hg38.pancan.sort.filt150.interval.bed merge.F5.hg38.pancan.sort.filt150.interval.bed merge.F5.hg38.pancan.sort.filt150.interval.bed > merge_3.bed

cat merge_3.bed | awk 'BEGIN{OFS="\t";} {print $1,$2,$3;}' > merge_3.tab.bed

bedtools shuffle -i merge_3.tab.bed \
 -g hg38.chrom.sizes \
 -excl repeat_hg38 \
 -seed 4242 \
 -noOverlapping > random_hg38.bed

sort -k1,1 -k2,2n random_hg38.bed > random_hg38.sort.bed

bedtools intersect -wa -wb \
 -a random_hg38.sort.bed \
 -b exon_hg38.bed hg38-blacklist.v2_sort.bed merge.F5.hg38.pancan.sort.filt150.interval.bed \
 -sorted \
 -v > filt_random_hg38.sort.bed

cat filt_random_hg38.sort.bed merge.F5.hg38.pancan.sort.filt150.bed | \
 awk 'BEGIN{OFS="\t";} {print $1,$2,$3;}' | sort -k1,1 -k2,2n > merge_filt_random_F5_pancan.hg38.sort.bed

