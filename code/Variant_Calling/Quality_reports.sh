#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --job-name=Quali_map
#SBATCH --output=Quali_map_%j.out
#SBATCH --error=Quali_map_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=120:00:00  # Adjusted walltime to 24 hours
#SBATCH --ntasks=48  # Using all available cores on a single node
#SBATCH --nodes=1  # Use a single node
#SBATCH --partition=atlas  # Selecting the atlas partition
#SBATCH --qos=normal  # Using the normal QoS
# Source the environment variables
source VC_project_env.sh

# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $VC_code"
echo "Data Directory: $VC_data"
echo "Results Directory: $VC_results"

# Define the main output directory where all job-specific directories are located
output_dir=${VC_results}/Quality_Control
mkdir -p ${output_dir}
cd ${output_dir}
mkdir -p ${output_dir}/original_mapping
mkdir -p ${output_dir}/dedup_mapping
mkdir -p ${output_dir}/pre_trim/
mkdir -p ${output_dir}/post_trim/
# Define BAM directories
bam_dir=${VC_results}/2.Alignment/
dedup_bam_dir=${VC_results}/3.Cleaning/
pre_fastq=${VC_results}/1.Trim_QC/pre_trim
post_fastq=${VC_results}/1.Trim_QC/post_trim
# Set Java options
export JAVA_OPTS="-Xms1G -Xmx8G"

# Load required modules
ml python
ml miniconda3


for report in ${bam_dir}/*.sorted_aligned_reads_stats/; do
    { 
    name=$(basename $report .sorted_aligned_reads_stats )
    echo $name
    cp -r ${report}/ ${output_dir}/original_mapping/   
    } 
done


for report in ${dedup_bam_dir}/*.marked_duplicates_stats/; do
    { 
    name=$(basename $report .marked_duplicates_stats )
    echo $name
    cp -r ${report}/ ${output_dir}/dedup_mapping/  
    } 
done



for report in ${pre_fastq}/*fastqc.zip; do
    { 
    cp -r ${report} ${output_dir}/pre_trim/   
    } 
done

for report in ${post_fastq}/*fastqc.zip; do
    { 
    cp -r ${report} ${output_dir}/post_trim/   
    } 
done


# Activate conda environment for Qualimap
eval "$(conda shell.bash hook)"
# Activate conda environment for MultiQC
conda activate ${my_softwares}/multiqc/env

# Run MultiQC
cd ${output_dir}
multiqc ./

# Deactivate conda environment
conda deactivate