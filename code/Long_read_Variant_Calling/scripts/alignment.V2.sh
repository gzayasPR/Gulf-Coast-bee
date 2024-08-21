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
ml minimap2
ml samtools
# Define variables
name=$1
reads=$2
ref_genome=$3
output_dir=$4

# Echo the inputs
echo "Inputs provided:"
echo "Sample Name: ${name}"
echo "Reads Directory: ${reads}"
echo "Reference Genome: ${ref_genome}"
echo "Output Directory: ${output_dir}"

#minimap2 -ax map-hifi ${ref_genome} ${reads}  > ${output_dir}/${name}.sam

echo "Start SAM to BAM"
samtools view -@ 48 -Sb  ${output_dir}/${name}.sam | samtools sort -@ 48 -o ${output_dir}/${name}.bam
samtools index -@ 48 ${output_dir}/${name}.bam
echo "End SAM to BAM"

# Generate alignment statistics
echo "Start qualimap"
ml miniconda3 
eval "$(conda shell.bash hook)"
source ~/.bashrc
# Check if trim_galore is available
conda activate ${my_softwares}/qualimap/env
export JAVA_OPTS="-Xms1G -Xmx8G"
qualimap bamqc -bam ${output_dir}/${name}.bam
echo "End qualimap"