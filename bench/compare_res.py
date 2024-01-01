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
        barcode = result_row['name'].split('_')[0]  # Assuming name is in the format BarcodeXX_96_test1
        barcode_number = int(barcode.replace('Barcode', ''))
        
        # Adjust barcode number for comparison if necessary
        if barcode_number >= 79:
            barcode_number += 1  # Adjusting for offset
        
        # Find the corresponding row in the baseline dataframe
        baseline_row = baseline_df[baseline_df['Nr'] == barcode_number]
        
        # If baseline_row is empty, no match was found in the baseline
        if baseline_row.empty:
            print(f"No baseline entry found for {barcode}")
            continue
              
        # Initialize match counters
        exact_matches = 0
        similar_matches = 0
        mismatches = 0
        not_classified = 0
        
        # Iterate over each locus
        for locus in ['A', 'B', 'C', 'DRB1', 'DQB1', 'DQA1', 'DPB1']:
            # Get results alleles
            result_alleles = result_row[[f'{locus}', f'{locus}.1']].fillna('').values  # .1 for the second allele in the result
            # Get baseline alleles
            baseline_alleles = list(baseline_row[[f'HLA-{locus}-1', f'HLA-{locus}-2']].fillna('').values.flatten())  # Assuming the baseline format
            
            # Clean up baseline alleles
            for num in range(len(baseline_alleles)): 
                if len(baseline_alleles[num]) > 3:
                    allele = re.sub('.+\s+','', baseline_alleles[num])
                    baseline_alleles[num]= allele.replace(f'{locus}*', '')
            
            # Perform comparison
            for result_allele in result_alleles:
                if barcode == 'Barcode30':
                    continue
                if locus == 'DQA1':
                    # For DQA1, check only the group number for an exact match
                    matched_alleles = [allele for allele in baseline_alleles if result_allele.split(':')[0] == allele]
                    if matched_alleles:
                        exact_matches += 1
                        baseline_alleles.remove(matched_alleles[0])
                    else:
                        mismatches += 1
                else:
                    # Exact match check for loci other than DQA1
                    if result_allele in baseline_alleles:
                        exact_matches += 1
                        baseline_alleles.remove(result_allele)
                    # Similar match check for loci other than DQA1
                    elif any(result_allele.split(':')[0] == baseline_allele.split(':')[0] for baseline_allele in baseline_alleles):
                        similar_matches += 1
                        for allele in baseline_alleles:
                            if result_allele.split(':')[0] == allele.split(':')[0]:
                                baseline_alleles.remove(allele)
                                break
                    # Mismatches and not classified
                    else:
                        if 'G' in result_allele:
                            mismatches += 1
                        else:
                            not_classified += 1
        
        # Store the results for this barcode
        comparison_results.append({
            'Barcode': barcode,
            'Exact Matches': exact_matches,
            'Similar Matches': similar_matches,
            'Mismatches': mismatches,
            'Not Classified': not_classified
        })
    # Convert the comparison results to a dataframe
    comparison_df = pd.DataFrame(comparison_results)
    # Calculate summary statistics and add them as a new row
    summary_row = {
        'Barcode': 'Summary',
        'Exact Matches': comparison_df['Exact Matches'].sum(),
        'Similar Matches': comparison_df['Similar Matches'].sum(),
        'Mismatches': comparison_df['Mismatches'].sum(),
        'Not Classified': comparison_df['Not Classified'].sum()
    }

    # Append the summary row to the dataframe
    summary_df = pd.DataFrame([summary_row], index=[comparison_df.index[-1] + 1])
    comparison_df = pd.concat([comparison_df, summary_df])

    # Save the comparison dataframe to a CSV file
    comparison_df.to_csv(output_file, index=False)

if __name__ == '__main__':
    # Taking command line arguments for the result file, baseline file, and output file
    result_file_path = sys.argv[1]
    baseline_file_path = sys.argv[2]
    output_file_path = sys.argv[3]
    
    # Perform comparison
    compare_alleles(result_file_path, baseline_file_path, output_file_path)
