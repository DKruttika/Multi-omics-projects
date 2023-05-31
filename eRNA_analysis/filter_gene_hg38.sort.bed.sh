#!/bin/bash

awk '{ if ($13 == "\"protein_coding\";") { print $0} }' gene_hg38.sort.bed > gene_hg38.sort.protein.coding.bed
