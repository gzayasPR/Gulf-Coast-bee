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
ml picard

# Define variables
name=$1
trimmed_dir=$2
R1=$3
R2=$4
ref_genome=$5
output_dir=$6
my_softwares=$7
# Unzip fastq.gz files
# Run bwa mem with 48 threads
echo "Start BWA MEM"
bwa mem -t 48 ${ref_genome} ${R1} ${R2} > ${output_dir}/${name}.aligned_reads.sam
echo "End BWA MEM"
# Convert SAM to BAM, sort and index
echo "Start SAM to BAM"
samtools view -@ 48 -Sb ${output_dir}/${name}.aligned_reads.sam | samtools sort -@ 48 -o ${output_dir}/${name}.sorted_aligned_reads.bam
samtools index -@ 48 ${output_dir}/${name}.sorted_aligned_reads.bam
echo "End SAM to BAM"

# Generate alignment statistics
echo "Start qualimap"
ml miniconda3 
eval "$(conda shell.bash hook)"
source ~/.bashrc
# Check if trim_galore is available
conda activate ${my_softwares}/qualimap/env
export JAVA_OPTS="-Xms1G -Xmx8G"
qualimap bamqc -bam ${output_dir}/${name}.sorted_aligned_reads.bam 
echo "End qualimap"