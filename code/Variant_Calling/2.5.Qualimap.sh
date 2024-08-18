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
source VC_project_env.sh
# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $VC_code"
echo "Data Directory: $VC_data"
echo "Results Directory: $VC_results"
# Define the main output directory where all job-specific directories are located
output_dir=$VC_results/Quality_Control/mapping/
cd ${output_dir}
bam_dir=
dedup_bam_dir=
export JAVA_OPTS="-Xms1G -Xmx8G"
ml miniconda3 
eval "$(conda shell.bash hook)"
conda activate ${my_softwares}/qualimap/env
qualimap bamqc -bam *.marked_duplicates.bam
conda deactivate
conda activate ${my_softwares}/multiqc/env
multiqc ./
conda deactivate

    