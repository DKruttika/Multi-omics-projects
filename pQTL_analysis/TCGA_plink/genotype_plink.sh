#!/bin/bash
#
#$ -cwd
#$ -pe smp 16
#$ -l mem_free=50G
#
#$ -S /bin/bash
#$ -e /common/bermanblab/data/private_data/POC_ROC/POC_ROC_pQTL/temp_plink_output/plink_tcga
#$ -o /common/bermanblab/data/private_data/POC_ROC/POC_ROC_pQTL/temp_plink_output/plink_tcga

source /home/dabkek/miniconda3/etc/profile.d/conda.sh
conda activate /home/dabkek/miniconda3/envs/pQTL

plink --bfile tcga_genotypes_cptac --pheno prot_data_tcga_pqtl_60_91721.txt --all-pheno --linear \
--allow-extra-chr --no-parents --allow-no-sex
