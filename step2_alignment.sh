#!/bin/bash

# Define input files and reference genome
FASTQ1="clipped_FC560466_S215_R1_001.fastq.gz"
FASTQ2="clipped_FC560466_S215_R2_001.fastq.gz"
REFERENCE="GAR104.fa"
OUTPUT_BAM="output.bam"

# Index the reference genome
bwa index $REFERENCE

# Perform alignment with BWA
bwa mem -R "@RG\tID:sample\tSM:sample\tLB:lib\tPL:illumina" \
-A 1 -B 4 -O 6,6 -E 1,1 -L 5,5 -U 17 \
$REFERENCE $FASTQ1 $FASTQ2 > $OUTPUT_BAM
