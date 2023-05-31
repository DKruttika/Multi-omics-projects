#!/bin/bash

bedtools intersect -wao -f 0.1\
 -a pancan_prognostic_eRNA_tab_sort.bed \
 -b sign_enh_calc_25k_64k_50923_tab_sort.bed \
 -sorted > common_pancan_POCROC_sort.bed
