#!/usr/bin/env python

import argparse
import csv
from collections import defaultdict

def parse_g_group_file(file_path):
    g_group_mapping = {}
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith("#") and line.strip():
                parts = line.strip().split(';')
                if parts[2] != '':  # ensuring that the G group exists
                    alleles = parts[1].strip().split('/')
                    g_group = parts[0].strip() + parts[2].strip()
                    for allele in alleles:
                        allele = parts[0].strip() + allele
                        g_group_mapping[allele] = g_group
    return g_group_mapping

def translate_to_g_groups(allele_results, g_group_mapping):
    translated_results = defaultdict()
    print(g_group_mapping.get('A*02:67'))
    for name, allele in allele_results:
        # print(allele)
        g_group = g_group_mapping.get(allele, allele)  # Default to original allele if no G group found
        translated_results[name] = g_group
    # print(translated_results)
    return translated_results

def read_results_file(file_path):
    results = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("HLA"):
                line = line.split(' ')
                results.append((line[0], line[1]))  # Assuming the format is [name, allele]
    return results

def write_output_file(translated_results, output_file):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        for name, allele in translated_results.items():
            writer.writerow([name, allele])

def main(args):
    g_group_mapping = parse_g_group_file(args.g_group_file)
    allele_results = read_results_file(args.input_file)
    translated_results = translate_to_g_groups(allele_results, g_group_mapping)
    write_output_file(translated_results, args.output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate HLA alleles to G groups")
    parser.add_argument('-i', '--input_file', required=True, help='Input file with HLA alleles')
    parser.add_argument('-g', '--g_group_file', required=True, help='File containing G group definitions')
    parser.add_argument('-o', '--output_file', required=True, help='Output file to write translated alleles')
    
    args = parser.parse_args()
    # args.g_group_file = '/Users/klarakotkowska/thesis/projectNEXTFLOW/data/hla_nom_g.txt'
    # args.input_file = '/Users/klarakotkowska/thesis/projectNEXTFLOW/work/55/c640956df41c3a1e7b1ef0d43cb242/Barcode94_96_test1_aggregated.txt'
    main(args)
