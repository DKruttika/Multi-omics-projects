#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l mem_free=20G
#$ -o /common/dabkek/enhancer_kallisto/diff_enh_quant/scratch
#$ -e /common/dabkek/enhancer_kallisto/diff_enh_quant/scratch

module load bedtools

bedtools getfasta -fi /common/bermanblab/data/public_data/refgenome/hg38.fa \
 -bed gene_transcript_blcklist_f5_hg38_og_sort.bed -fo filt_og_f5_hg38.fa.out


bedtools getfasta -fi /common/bermanblab/data/public_data/refgenome/hg38.fa \ 
 -bed gene_transcript_blcklist_f5_hg38_6kb_sort.bed -fo filt_6kb_f5_hg38.fa.out
