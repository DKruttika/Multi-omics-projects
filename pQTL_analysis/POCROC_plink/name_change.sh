#!/bin/bash


INFILE="$1.assoc.linear"
OUTFILE="./temp_plink_output/$1_pvalue.txt"

awk 'NR == 1{print $9 "|" FILENAME;next;} {print $9}' ${INFILE} > ${OUTFILE}
