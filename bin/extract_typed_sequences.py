#!/usr/bin/env python

import sys

def parse_name(line):
    # Extracts the allele name part, assuming it's the second part in a space-separated string
    parts = line.split()
    return parts[1] if len(parts) > 1 else None

def read_names(names_file):
    with open(names_file, 'r') as file:
        return set(parse_name(name.strip()) for name in file)

def extract_sequences(fasta_file, names_set, output_file):
    with open(fasta_file, 'r') as fasta, open(output_file, 'w') as output:
        write_sequence = False
        for line in fasta:
            if line.startswith('>'):
                allele_name = parse_name(line[1:])  # Remove '>' and parse
                write_sequence = allele_name in names_set if allele_name else False
            if write_sequence:
                output.write(line.upper())

if __name__ == '__main__':
    names_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]

    names_set = read_names(names_file)
    extract_sequences(fasta_file, names_set, output_file)
