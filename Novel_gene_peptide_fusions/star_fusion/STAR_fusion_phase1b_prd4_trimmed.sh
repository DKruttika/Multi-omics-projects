#!/bin/bash
#
#$ -cwd
#$ -pe smp 10
#$ -l mem_free=95G,h_vmem=140G
#$ -S /bin/bash
#$ -e path_to_working_dir
#$ -o path_to_working_dir


set -xe

INDIR="/common/bermanblab/data/private_data/POC_ROC/coetzeesg_pocroc/phase1/fastq_refix/RNA2-MJ-POCROC-20180306"
INFILE1="${INDIR}/$(basename ${INDIR}_R1.fastq.gz)"
INFILE2="${INDIR}/$(basename ${INDIR}_R2.fastq.gz)"
OUTDIR="/common/bermanblab/data/private_data/POC_ROC/POC_ROC_RNA-Seq/Phase1_star_fusion/$(basename ${INDIR})"
CORES="10"

mkdir -p $OUTDIR

[[ -e "${INFILE1}" ]] && echo "infile1 ${INFILE1} exists"
[[ -e "${INFILE2}" ]] && echo "infile2 ${INFILE2} exists"


conda activate starfusion

STAR-Fusion \
	--genome_lib_dir /common/bermanblab/data/private_data/POC_ROC/coetzeesg_pocroc/reference/GRCh38_gencode_v29_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir \
	--left_fq ${INFILE1}\
	--right_fq ${INFILE2} \
	--FusionInspector validate \
	--examine_coding_effect \
	--denovo_reconstruct \
	--CPU ${CORES} \
	--output_dir ${OUTDIR}
