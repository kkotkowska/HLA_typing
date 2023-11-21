#!/usr/bin/env python

import csv
import argparse
from collections import defaultdict

def extract_locus_and_type(line):
    """ Extract gene type and number from the line. """
    parts = line.split(',')
    locus = parts[1].split('*')[0]
    type = parts[1].split('*')[1].split()[0]
    return locus, type

def process_file(input_file, output_file):
    data = defaultdict(list)

    # Read and process the input file
    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith("HLA"):
                locus, type = extract_locus_and_type(line)
                data[locus].append(type)

    # Sort the data
    sorted_data = sorted(data.items(), key=lambda x: (x[0], x[1]))

    # Write to the output CSV file
    with open(output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Locus", "Gene Type"])

        for locus, types in sorted_data:
            if len(types) == 1:  # Duplicate if only one record
                types *= 2
            for type in sorted(types):
                csvwriter.writerow([locus, type])

def main():
    parser = argparse.ArgumentParser(description="Process HLA data into CSV format.")
    parser.add_argument('input_file', type=str, help="Path to the input text file.")
    parser.add_argument('output_file', type=str, help="Path to the output CSV file.")

    args = parser.parse_args()
    process_file(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
