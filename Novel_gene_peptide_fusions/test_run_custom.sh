#!/bin/bash
#
#$ -cwd
#$ -pe smp 6
#$ -l mem_free=10G
#$ -S /bin/bash



conda activate FusionPro
perl -I ${PWD}/bin ${PWD}/bin/FusionPro-C.pl -bin ${PWD}/bin -filter IP -cdna ${PWD}/Homo_sapiens.GRCh38.cdna.82.all.fa -gtf ${PWD}/exon_coord_for_custom.tsv -i ${PWD}/fusionPro_input1.txt -o ${PWD} -r ${PWD}/library
