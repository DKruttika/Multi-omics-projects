#!/bin/bash
#
#$ -cwd
#$ -pe smp 6
#$ -l mem_free=30G
#$ -S /bin/bash


#set -xe

#INDIR="/common/bermanblab/data/private_data/POC_ROC/coetzeesg_pocroc/phase2/STAR_trimmed/FT-SA13637_S1_merged"
#INFILE1="${INDIR}/$(basename ${INDIR}_Aligned.out.sorted.bam)"
#OUTDIR="/common/bermanblab/data/private_data/POC_ROC/POC_ROC_RNA-Seq/arriba_fusions/phase2/$(basename ${INDIR})"

#mkdir -p $OUTDIR

#[[ -e "${INFILE1}" ]] && echo "infile1 ${INFILE1} exists"


#INFILE1=$1

conda activate arriba

ARRIBA_FILES=/home/dabkek/miniconda3/envs/arriba/var/lib/arriba

/common/bermanblab/data/private_data/POC_ROC/POC_ROC_RNA-Seq/arriba_fusions/run_arriba_on_prealigned_bam.sh /common/bermanblab/data/private_data/POC_ROC/POC_ROC_RNA-Seq/arriba_fusions/STAR_index_hg38_GENCODE28/ \
/common/bermanblab/data/private_data/POC_ROC/POC_ROC_RNA-Seq/arriba_fusions/GENCODE28.gtf \
/common/bermanblab/data/private_data/POC_ROC/POC_ROC_RNA-Seq/arriba_fusions/hg38.fa \
$ARRIBA_FILES/blacklist_hg38_GRCh38_v2.2.1.tsv.gz \
$ARRIBA_FILES/known_fusions_hg38_GRCh38_v2.2.1.tsv.gz \
$ARRIBA_FILES/protein_domains_hg38_GRCh38_v2.2.1.gff3 \
8 \
/common/bermanblab/data/private_data/POC_ROC/coetzeesg_pocroc/phase1/rsem_trimmed_repair/.genome.sorted.bam

