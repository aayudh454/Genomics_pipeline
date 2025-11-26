#!/bin/bash

# Define input BAM and output sorted BAM
INPUT_BAM="output.bam"                     # Input BAM file (aligned reads)
SORTED_BAM="sorted_output.bam"             # Output BAM file (sorted reads)

# Sort and index using samtools
samtools sort $INPUT_BAM -o $SORTED_BAM    # Sort the BAM file and save the sorted output to $SORTED_BAM
samtools index $SORTED_BAM                 # Index the sorted BAM file for efficient access

