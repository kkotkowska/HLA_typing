import pandas as pd
import sys
import re

def compare_alleles(result_file, baseline_file, output_file):
    # Load the results file
    results_df = pd.read_csv(result_file)
    
    # Load the baseline Excel file
    baseline_df = pd.read_excel(baseline_file, engine='openpyxl')

    # Initialize a list to hold comparison results
    comparison_results = []

    # Iterate over each row in the result dataframe
    for index, result_row in results_df.iterrows():
        barcode = result_row['name'].split('_')[0]
        barcode_number = int(barcode.replace('Barcode', ''))

        if barcode_number >= 79:
            barcode_number += 1
        if barcode_number == 30:
            continue
        
        # Find the corresponding row in the baseline dataframe
        baseline_row = baseline_df[baseline_df['Nr'] == barcode_number]

        if baseline_row.empty:
            print(f"No baseline entry found for {barcode}")
            continue
        
        # Initialize counters for each comparison case
        both_correct = 0
        both_incorrect = 0
        heterozygous_instead_of_homozygous = 0
        homozygous_instead_of_heterozygous = 0
        one_incorrect = 0

        # Iterate over each locus
        for locus in ['A', 'B', 'C', 'DRB1', 'DQB1', 'DQA1', 'DPB1']:
            # Get results alleles
            result_alleles = result_row[[f'{locus}', f'{locus}.1']].fillna('').values
            
            # Get baseline alleles
            baseline_alleles = baseline_row[[f'HLA-{locus}-1', f'HLA-{locus}-2']].fillna('').values.flatten()
            baseline_alleles = [allele.split()[-1] for allele in baseline_alleles if allele]
            
            # Clean up baseline alleles
            for num in range(len(baseline_alleles)): 
                if len(baseline_alleles[num]) > 3:
                    allele = re.sub('.+\s+','', baseline_alleles[num])
                    baseline_alleles[num]= allele.replace(f'{locus}*', '')
            
            # Check for comparison cases
            # print(result_alleles)
            # print(baseline_alleles)
            # Check for homozygous, heterozygous, and incorrect calls
            if locus == 'DQA1':
                # Special handling for DQA1
                result_allele_numbers = [allele.split(':')[0] for allele in result_alleles]
                baseline_allele_numbers = [allele.split(':')[0] for allele in baseline_alleles]
                if len(set(result_allele_numbers)) == 1 and len(set(baseline_allele_numbers)) == 1 and result_allele_numbers == baseline_allele_numbers:
                    # Homozygous
                    both_correct += 1
                elif len(set(result_allele_numbers)) == 2 and len(set(baseline_allele_numbers)) == 2 and result_allele_numbers == baseline_allele_numbers:
                    # Heterozygous
                    both_correct += 1
                elif len(set(result_allele_numbers)) == 1 and len(set(baseline_allele_numbers)) == 2:
                    homozygous_instead_of_heterozygous += 1
                elif len(set(result_allele_numbers)) == 2 and len(set(baseline_allele_numbers)) == 1:
                    heterozygous_instead_of_homozygous += 1
                elif not any(num in baseline_allele_numbers for num in result_allele_numbers):
                    both_incorrect += 1
                else:
                    one_incorrect += 1
                
            else:
                if len(set(result_alleles)) == 1 and len(set(baseline_alleles)) == 1 and result_alleles[0] in baseline_alleles:  # Called homozygous correctly
                    both_correct += 1
                elif len(set(result_alleles)) == 2 and all(allele in baseline_alleles for allele in result_alleles):  # Called heterozygous correctly
                    both_correct += 1
                elif len(set(result_alleles)) == 1 and not any(allele in baseline_alleles for allele in result_alleles):  # Called homozygous incorrectly
                    both_incorrect += 1
                elif len(set(result_alleles)) == 2 and not any(allele in baseline_alleles for allele in result_alleles):  # Called heterozygous incorrectly
                    both_incorrect += 1
                elif len(set(result_alleles)) == 1 and len(set(baseline_alleles)) == 2:  # Called homozygous when it is heterozygous
                    homozygous_instead_of_heterozygous += 1
                elif len(set(result_alleles)) == 2 and len(set(baseline_alleles)) == 1:  # Called heterozygous when it is homozygous
                    heterozygous_instead_of_homozygous += 1
                else:
                    one_incorrect += 1

        # Store the results for this barcode
        comparison_results.append({
            'Barcode': barcode,
            'Both Alleles Correct': both_correct,
            'Both Alleles Incorrect': both_incorrect,
            'Heterozygous Instead of Homozygous': heterozygous_instead_of_homozygous,
            'Homozygous Instead of Heterozygous': homozygous_instead_of_heterozygous,
            'Zygote Correct - One Allele Incorrect': one_incorrect
        })

    # Convert the comparison results to a dataframe
    comparison_df = pd.DataFrame(comparison_results)

    # Save the comparison dataframe to a CSV file
    comparison_df.to_csv(output_file, index=False)

if __name__ == '__main__':
    result_file_path = sys.argv[1]
    baseline_file_path = sys.argv[2]
    output_file_path = sys.argv[3]
    
    compare_alleles(result_file_path, baseline_file_path, output_file_path)

