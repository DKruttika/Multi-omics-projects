#!/bin/bash

bedtools intersect -wao -f 0.1\
 -a GeneHancer_AnnotSV_elements_v5.15.tab.sort.bed \
 -b sign_enh_275.tab.sort.bed \
 -sorted > common_genehancer_sign_erna.sort.bed
