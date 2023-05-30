#!/bin/bash
#
#$ -cwd
#$ -pe smp 6
#$ -l mem_free=50G
#
#$ -S /bin/bash
#$ -e /common/bermanblab/data/private_data/POC_ROC/POC_ROC_pQTL/scratch
#$ -o /common/bermanblab/data/private_data/POC_ROC/POC_ROC_pQTL/scratch

source /home/dabkek/miniconda3/etc/profile.d/conda.sh
conda activate /home/dabkek/miniconda3/envs/pQTL

plink --vcf test_mac_6_MJ.vcf --pheno prot_data_noimp_filt_20_pqtl_ph1_ph2_60821.txt --all-pheno --linear \
--allow-extra-chr --no-parents --allow-no-sex --maf 0.1
