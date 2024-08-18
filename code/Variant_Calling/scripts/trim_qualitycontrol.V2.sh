#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --job-name="runQC_$1"
#SBATCH --output=runQC_$1_%j.out
#SBATCH --error=runQC_$1_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=24:00:00  # Adjusted walltime to 24 hours
#SBATCH --ntasks=48  # Using all available cores on a single node
#SBATCH --nodes=1  # Use a single node
#SBATCH --partition=atlas  # Selecting the atlas partition
#SBATCH --qos=normal  # Using the normal QoS

name=$1
R1=$2
R2=$3
output_dir=$4
my_softwares=$5
ml fastqc
fastqc -t 48  ${R1} ${R2} -o ${output_dir}/pre_trim


input_file1=${R1}
input_file2=${R2}
# Define the input and output files
output_paired_file1=${output_dir}/post_trim/${name}_R1.paired.fastq
output_unpaired_file1=${output_dir}/post_trim/${name}_R1.unpaired.fastq
output_paired_file2=${output_dir}/post_trim/${name}_R2.paired.fastq
output_unpaired_file2=${output_dir}/post_trim/${name}_R2.unpaired.fastq

# Run Trime Galore
ml miniconda3 
eval "$(conda shell.bash hook)"
source ~/.bashrc
# Check if trim_galore is available
conda activate ${my_softwares}/trim_galore/env

trim_galore -j 4 --quality 30 --fastqc --length 80  -clip_R1 15 -clip_R2 15 --three_prime_clip_R1 3 \
            --three_prime_clip_R2 3 -o ${output_dir}/post_trim/ --paired ${input_file1} ${input_file2} --basename ${name} --dont_gzip

conda deactivate
