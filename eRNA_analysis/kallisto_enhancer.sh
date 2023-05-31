#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=24G
#$ -pe smp 8
#$ -S /bin/bash
#$ -o /common/dabkek/enhancer_kallisto/diff_enh_quant/scratch
#$ -e /common/dabkek/enhancer_kallisto/diff_enh_quant/scratch

set -xe

module load kallisto

#first input should be reads_1.fastq.gz- entire file name


INFILE1="/common/bermanblab/data/private_data/POC_ROC/coetzeesg_pocroc/phase1/fastq_refix/$1/$1_R1.fastq.gz"
INFILE2="/common/bermanblab/data/private_data/POC_ROC/coetzeesg_pocroc/phase1/fastq_refix/$1/$1_R2.fastq.gz"
OUTDIR="$1"
REF="merge_filt_random_F5_pancan.hg38.fa.out.idx"

#OUTDIR="/common/dabkek/test_kallisto/test"
#REF="/common/dabkek/test_kallisto/transcripts.fasta.gz"
#CORES="8"


[[ -e "${INFILE1}" ]] && echo "infile1 ${INFILE1} exists"
[[ -e "${INFILE2}" ]] && echo "infile2 ${INFILE2} exists"

mkdir -p $OUTDIR

#a="/hpc/apps/kallisto/0.44/kallisto"

#kallisto index -i 6_transcripts.idx /common/dabkek/test_kallisto/test/transcripts.fasta.gz

kallisto quant \
	--index="$(readlink -f ${REF})" \
	--output-dir="$(readlink -f ${OUTDIR})" \
	--rf-stranded \
	--bias \
	--bootstrap-samples=100 \
	--seed=4242 \
	"${INFILE1}" "${INFILE2}"
