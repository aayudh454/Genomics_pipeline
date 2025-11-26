#!/bin/bash
  
# Define input FASTQ files
FASTQ1="FC560466_S215_R1_001.fastq.gz"
FASTQ2="FC560466_S215_R2_001.fastq.gz"

# Run fastp with the specified parameters
fastp --detect_adapter_for_pe \
-W 4 \
-M 20 \
--cut_tail \
--length_required 15 \
-i $FASTQ1 \
-I $FASTQ2 \
-o clipped_$FASTQ1 \
-O clipped_$FASTQ2
