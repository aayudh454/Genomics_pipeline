#!/usr/bin/env python3 
import pandas as pd

def split_biallelic_rows(df):
    """
    Split rows with multiple alternate alleles into separate rows,
    each with its own alternate allele frequency.
    
    Parameters:
    df (DataFrame): Input dataframe.
    
    Returns:
    DataFrame: Modified dataframe with split alternate alleles.
    """
    rows = []
    for _, row in df.iterrows():
        alts = row['ALT'].split(',')
        ad_values = list(map(int, row['sample.AD'].split(',')))

        total_reads = sum(ad_values)
        for i, alt in enumerate(alts, start=1):
            new_row = row.copy()
            new_row['ALT'] = alt
            new_row['alt_freqHQ'] = ad_values[i] / total_reads if total_reads != 0 else 0
            new_row['var_type'], bp_diff = determine_variant_type(row['REF'], alt)
            new_row['bp_differences'] = bp_diff
            new_row['Frameshift'] = is_frameshift(new_row['var_type'], bp_diff)
            rows.append(new_row)

    return pd.DataFrame(rows)


def determine_variant_type(ref, alt):
    """
    Determine the type of variant and calculate base pair differences.

    Parameters:
    ref (str): The reference allele.
    alt (str): The alternate allele.

    Returns:
    tuple: A tuple containing:
        - var_type (str): 'SUB', 'INS', 'DEL', or 'MIXED'
        - bp_differences (int): Base pair difference for the ALT allele.
    """
    if alt == '*':
        return 'DEL', -len(ref)
    elif len(ref) == len(alt):
        return 'SUB', 0
    elif len(ref) > len(alt):
        return 'DEL', len(alt) - len(ref)  # Negative value
    else:
        return 'INS', len(alt) - len(ref)  # Positive value


def is_frameshift(var_type, bp_difference):
    """
    Determine if the variant causes a frameshift.

    Parameters:
    var_type (str): The type of variant ('INS', 'DEL', 'SUB', 'MIXED').
    bp_difference (int): Base pair difference.

    Returns:
    str: 'TRUE' if bp_difference is not a multiple of 3, 'FALSE' otherwise.
    """
    if var_type in ['INS', 'DEL']:
        return 'TRUE' if abs(bp_difference) % 3 != 0 else 'FALSE'
    return 'FALSE'


def process_vcf_table(input_file, output_file):
    """
    Process the input VCF table (tab-delimited) and split multi-allelic rows.

    The following columns are added/modified:
    - alt_freqHQ: Alternate allele frequency calculated from sample.AD.
    - altHQ_filt: 'TRUE' if alt_freqHQ >= 0.85, 'FALSE' otherwise.
    - var_type: The type of variant (SUB, INS, DEL, MIXED) derived from REF and ALT.
    - bp_differences: Base pair differences between REF and ALT.
    - Frameshift: 'TRUE' if the variant causes a frameshift, 'FALSE' otherwise.
    - Filter1: 'TRUE' if both altHQ_filt and Frameshift are 'TRUE', else 'FALSE'.

    Parameters:
    input_file (str): Path to the input tab-delimited file containing VCF data.
    output_file (str): Path to the output CSV file with additional columns.

    Returns:
    None
    """
    # Read the input file (tab-delimited)
    df = pd.read_csv(input_file, sep='\t')

    # Split rows with bi-allelic variants into separate rows
    df = split_biallelic_rows(df)

    # Calculate altHQ_filt
    df['altHQ_filt'] = df['alt_freqHQ'].apply(lambda x: 'TRUE' if x >= 0.85 else 'FALSE')

    # Calculate Filter1
    df['Filter1'] = df.apply(lambda row: 'TRUE' if row['altHQ_filt'] == 'TRUE' and row['Frameshift'] == 'TRUE' else 'FALSE', axis=1)

    # Save the updated dataframe to a final CSV file
    df.to_csv(output_file, index=False)


if __name__ == "__main__":
    # Define input and output files
    input_file = 'output.table'
    output_file = 'final_output.csv'

    # Process the VCF table and create the final output
    process_vcf_table(input_file, output_file)

