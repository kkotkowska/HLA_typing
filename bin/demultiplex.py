#!/usr/bin/env python

from Bio import SeqIO
import csv
import argparse
import os
import gzip
from Bio import Align

def smith_waterman_alignment(primer, sequence):
    """
    Perform Smith-Waterman alignment between primer and sequence using Biopython's pairwise aligner.
    Returns the alignment score.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"  # Local alignment
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1

    alignments = aligner.align(primer, sequence)
    return alignments.score

def reverse_slicing(s):
    return s[::-1]


# Define a function to demultiplex the FASTQ file based on primer sequences
def demultiplex_fastq(input_fastq, primer_csv, output_directory):
    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)
    # Create a dictionary to store output file handles for each primer prefix
    prefix_handles = {}
    # Create a dictionary to store primer sequences for each primer prefix
    primer_prefixes = {}

    # Read the CSV file containing primer information
    with open(primer_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=';')
        for row in reader:
            primer_name = row['PrimerName']
            primer_prefix = primer_name.split('-')[0]  # Extract the prefix
            primer_sequence = row['PrimerSequence'].upper()

            # Open an output file for each unique prefix
            if primer_prefix not in prefix_handles:
                output_file_path = os.path.join(output_directory, f"{primer_prefix}.fastq")
                prefix_handles[primer_prefix] = open(output_file_path, 'w')
            
            # Add the primer sequence to the primer_prefixes dictionary
            if primer_prefix not in primer_prefixes:
                primer_prefixes[primer_prefix] = []
            primer_prefixes[primer_prefix].append(primer_sequence)

    # Define a threshold for the acceptable alignment score
    alignment_threshold = len(primer_sequence) * 2 * 0.7

    # Iterate through the FASTQ file and write reads to the corresponding output files
    with gzip.open(input_fastq, 'rt') as fastq_file:
        for record in SeqIO.parse(fastq_file, 'fastq'):
            found_primer = False
            for prefix, primers in primer_prefixes.items():
                for primer_sequence in primers:
                    reversed_primer_sequence = reverse_slicing(primer_sequence)
                    if (smith_waterman_alignment(primer_sequence, str(record.seq)[:len(primer_sequence)+10]) >= alignment_threshold or
                            smith_waterman_alignment(reversed_primer_sequence, str(record.seq)[-len(reversed_primer_sequence)-10:]) >= alignment_threshold):
                        prefix_handles[prefix].write(record.format("fastq"))
                        found_primer = True
                        break
                if found_primer:
                    break

            # If no matching primer was found, write the read to an "unclassified" file
            if not found_primer:
                unclassified_file.write(record.format("fastq"))

    # Close all output files
    for output_file in prefix_handles.values():
        output_file.close()



if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Demultiplex FASTQ data based on primer sequences.")
    parser.add_argument("input_fastq", help="Input FASTQ file")
    parser.add_argument("primer_csv", help="CSV file containing primer information")
    parser.add_argument("output_directory", help="Name of the output directory")
    # parser.add_argument("output_dir", help="Output directory for demultiplexed files")
    args = parser.parse_args()

    # unclassified_file_path = os.path.join(args.output_dir, "unclassified.fastq")
    unclassified_file_path = "unclassified.fastq"
    unclassified_file = open(unclassified_file_path, 'w')  # Output file for unclassified reads

    # Demultiplex the FASTQ file
    # demultiplex_fastq(args.input_fastq, args.primer_csv, args.output_dir)
    demultiplex_fastq(args.input_fastq, args.primer_csv, args.output_directory)


    unclassified_file.close()
