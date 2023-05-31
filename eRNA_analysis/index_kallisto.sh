#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#$ -l mem_free=20G
#$ -e /common/dabkek/enhancer_kallisto/diff_enh_quant/scratch
#$ -o /common/dabkek/enhancer_kallisto/diff_enh_quant/scratch

module load bedtools
module load kallisto

bedtools getfasta -fi /common/bermanblab/data/public_data/refgenome/hg38.fa \
 -bed merge_filt_random_F5_pancan.hg38.sort.bed -fo merge_filt_random_F5_pancan.hg38.fa.out

kallisto index -i merge_filt_random_F5_pancan.hg38.fa.out.idx merge_filt_random_F5_pancan.hg38.fa.out
