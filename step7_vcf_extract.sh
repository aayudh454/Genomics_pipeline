#!/bin/bash
  
# Extract fields from the VCF to a table, DP from INFO field
gatk VariantsToTable \
    -V output.vcf.gz \
    -F CHROM -F POS -F REF -F ALT -F QUAL -F DP \
    -GF GT -GF AD \
    -O output.table
