#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --job-name="runAlign_$1"
#SBATCH --output=runAlign_$1_%j.out
#SBATCH --error=runAlign_$1_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=120:00:00  # Adjusted walltime to 24 hours
#SBATCH --ntasks=48  # Using all available cores on a single node
#SBATCH --nodes=1  # Use a single node
#SBATCH --partition=atlas  # Selecting the atlas partition
#SBATCH --qos=normal  # Using the normal QoS

# Load necessary modules
ml bwa
ml samtools
ml qualimap


# Define variables
name=$1
trimmed_dir=$2
ref_genome=$3
output_dir=$4

mkdir -p ${output_dir}/${name}_quality.report
# Run bwa mem with 48 threads
echo "Start BWA MEM"
#bwa mem -t 48 ${ref_genome} ${trimmed_dir}/${name}_R1.paired.fastq ${trimmed_dir}/${name}_R2.paired.fastq > ${output_dir}/${name}.aligned_reads.sam
echo "End BWA MEM"
# Convert SAM to BAM, sort and index
echo "Start SAM to BAM"
#samtools view -@ 48 -Sb ${output_dir}/${name}.aligned_reads.sam | samtools sort -@ 48 -o ${output_dir}/${name}.sorted_aligned_reads.bam
#samtools index -@ 48 ${output_dir}/${name}.sorted_aligned_reads.bam
echo "End SAM to BAM"

# Generate alignment statistics
echo "Start flagstat"
samtools flagstat ${output_dir}/${name}.sorted_aligned_reads.bam > ${output_dir}/${name}_quality.report/${name}.alignment_stats.txt
echo "end flagstat"
# Calculate depth of coverage
#bedtools genomecov -ibam ${output_dir}/${name}.sorted_aligned_reads.bam -d > #${output_dir}/${name}.genome_coverage.txt
#awk '{sum+=$3; count++} END {print "Average coverage depth:", sum/count}' ${output_dir}/${name}.genome_coverage.txt > ${output_dir}/${name}.coverage_summary.txt
echo "Start qualimap"
export JAVA_OPTS="-Xms1G -Xmx8G"
# Assess alignment quality using Qualimap
qualimap bamqc -bam ${output_dir}/${name}.sorted_aligned_reads.bam -outdir ${output_dir}/${name}_quality.report/ --java-mem-size=8G --nt 48
echo "End qualimap"
