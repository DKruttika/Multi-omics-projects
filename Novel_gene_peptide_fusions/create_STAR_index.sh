#!/bin/bash
#
#$ -cwd
#$ -pe smp 6
#$ -l mem_free=30G
#$ -S /bin/bash


conda activate arriba

/home/dabkek/miniconda3/envs/arriba/var/lib/arriba/download_references.sh hg38+GENCODE28
