#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --job-name="runQC_$1"
#SBATCH --output=runQC_$1_%j.out
#SBATCH --error=runQC_$1_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=24:00:00  # Adjusted walltime to 24 hours
#SBATCH --ntasks=48  # Using all available cores on a single node
#SBATCH --nodes=2  # Use a single node
#SBATCH --partition=atlas  # Selecting the atlas partition
#SBATCH --qos=normal  # Using the normal QoS

# Echo the inputs
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

ml miniconda3 
eval "$(conda shell.bash hook)"
source ~/.bashrc
# Check if trim_galore is available
conda activate ${my_softwares}/fastqc/env
fastqc -nogroup -t 46  ${read} -o ${output_dir}

