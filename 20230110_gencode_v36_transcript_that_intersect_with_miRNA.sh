#!/bin/bash

# This script is to intersect transcript bed file and miRNA bed file and generate list about transcript that intersect with miRNA 

bedtools intersect -a miRBase_v21.bed -b UCSC_gencode_v36_comprehensive_transcript.bed -wa -wb -s > TCGA_hg38_transcript_intersect_with_miRNA.txt