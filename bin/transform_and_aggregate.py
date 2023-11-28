#!/usr/bin/env python

import csv
import sys
from collections import defaultdict

def read_csv(file_path):
    # Read a CSV file and return data in a desired format
    data = defaultdict(list)
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            locus = row['Locus']
            gene_type = row['Gene Type']
            data[locus].append(gene_type)
    return data

def extract_test_name(file_path):
    name = file_path.split('/')[-1].replace('_translated.csv', '')
    return name

def main(csv_files, output_file):
    aggregated_data = defaultdict(lambda: defaultdict(list))
    loci = ['A', 'B', 'C', 'DRB1', 'DQB1', 'DQA1', 'DPB1']  # Adjust this list based on your actual loci
    
    # Sort the csv_files list
    sorted_csv_files = sorted(csv_files, key=lambda x: extract_test_name(x))

    for file_path in sorted_csv_files:
        test_name = extract_test_name(file_path)  # Implement this function based on your filename pattern
        data = read_csv(file_path)
        for locus, gene_types in data.items():
            aggregated_data[test_name][locus] = gene_types

    # Now write aggregated_data to a CSV file in the desired format
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        # Write header
        header = ['name'] + [locus for locus in loci for _ in range(2)]  # Two columns for each locus
        writer.writerow(header)

        # Write rows
        for test_name, loci_data in aggregated_data.items():
            row = [test_name]
            for locus in loci:
                gene_types = loci_data.get(locus, [''] * 2)  # Ensure two entries for each locus
                row.extend(gene_types[:2])  # Limit to two entries per locus
            writer.writerow(row)

if __name__ == "__main__":
    csv_files = sys.argv[1:-1]
    output_file = sys.argv[-1]
    main(csv_files, output_file)
