#/bin/bash

FILE1="$1"
OUTPUT="./genotype_allele/$1"

awk '(FNR==NR){minor[FNR]=$1; het[FNR]=$2; major[FNR]=$3; next} (FNR!=NR && FNR==1){print; next} (FNR!=NR && $1==0){print minor[FNR-1]; next} (FNR!=NR && $1==1){print het[FNR-1]; next} (FNR!=NR && $1==2){print major[FNR-1]; next} (FNR!=NR && $1==-1){print "0"; next} (FNR!=NR){print "invalid"; next}' minor_major.txt ${FILE1} > ${OUTPUT}
