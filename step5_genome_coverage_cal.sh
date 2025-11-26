#!/bin/bash

# Define input BAM file
INPUT_BAM="sorted_output.bam"
COVERAGE_SUMMARY="coverage_summary.txt"

# Calculate genome coverage
bedtools genomecov -ibam $INPUT_BAM | grep "^genome" > $COVERAGE_SUMMARY
