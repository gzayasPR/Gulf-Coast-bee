#!/bin/bash

# Define the main output directory where all job-specific directories are located
output_dir=/90daydata/beenome100/hesperapis_oraria_genomics/Hesperapis_oraria_Pop_genomics/results/Variant_Calling/2.Alignment/
# Define the output summary file
summary_file="${output_dir}/alignment_quality_summary.txt"
rm $summary_file
# Initialize the summary file
echo -e "Sample\tMapped Reads\tProperly Paired Reads\tMean Coverage\tCoverage StdDev\tGC Content\tMapping Quality" > $summary_file

# Iterate over each job directory to extract metrics
for sample in ${output_dir}/*.aligned_reads.sam; do
    sample_name=$(basename $sample .aligned_reads.sam)
    echo $sample_name
    # Extract metrics from samtools flagstat
             flagstat_file="${output_dir}/${sample_name}_quality.report/${sample_name}.alignment_stats.txt"
    mapped_reads=$(grep -m 1 "mapped (" $flagstat_file | awk '{print $1}')
    mapped_reads=$mapped_reads$(grep -m 1 "mapped (" $flagstat_file | awk '{print $5}')")"
    pair_aligned=$(grep "properly paired" $flagstat_file | awk '{print $1}')$(grep "properly paired" $flagstat_file | awk '{print $5}')")"
    # Extract metrics from qualimap report
    qualimap_report="${output_dir}/${sample_name}_quality.report/genome_results.txt"
    mean_coverage=$(grep "mean coverageData =" $qualimap_report | awk '{print $4}')
    std_coverage=$(grep "std coverageData =" $qualimap_report | awk '{print $4}')
    gc_content=$(grep "GC percentage" $qualimap_report | awk '{print $4}')
    mapping_quality=$(grep "mean mapping quality" $qualimap_report | awk '{print $5}')
    # Append the metrics to the summary file
    echo -e "${sample_name}\t${mapped_reads}\t${pair_aligned}\t${mean_coverage}\t${std_coverage}\t${gc_content}\t${mapping_quality}" >> $summary_file
done

echo "Aggregation complete. Results written to: $summary_file"
