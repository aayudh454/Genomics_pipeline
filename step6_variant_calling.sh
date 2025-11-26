#!/bin/bash

# Define input files and reference genome
INPUT_BAM="sorted_output.bam"
REFERENCE="GAR104.fa"
VCF_OUTPUT="output.vcf.gz"

# Run HaplotypeCaller and compress output
gatk HaplotypeCaller \
    -R $REFERENCE \
    -I $INPUT_BAM \
    -O $VCF_OUTPUT \
    --native-pair-hmm-threads 32

# Index the VCF (compressed and forced overwrite)
tabix -f -p vcf $VCF_OUTPUT
