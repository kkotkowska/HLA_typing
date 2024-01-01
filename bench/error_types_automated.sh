#!/bin/bash

# Define an array of threshold values
threshold_values=(0.0 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
# threshold_values=(0.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1)

# Output directory prefix
# output_dir_prefix="res_s3s2"

# Loop through each threshold value
for threshold in "${threshold_values[@]}"
do
    # Define the output directory for each threshold
    # output_dir="${output_dir_prefix}_threshold_${threshold}"

    # Execute the Nextflow command with the current threshold value
    python benchmarking/error_types_analysis.py projectNEXTFLOW/results/s3s2_threshold_${threshold}.csv benchmarking/baseline.xlsx benchmarking/s3s2_${threshold}.csv
    # Optional: Add a delay or log message if needed
    # sleep 1
    # echo "Completed run for threshold $threshold"
done

# End of script
