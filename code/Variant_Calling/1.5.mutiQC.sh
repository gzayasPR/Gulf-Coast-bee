#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=MultiQC_%j.out
#SBATCH --error=MultiQC_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb

# Define project directories
source VC_project_env.sh

# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $VC_code"
echo "Data Directory: $VC_data"
echo "Results Directory: $VC_results"
source ~/.bashrc
multi_qc=${my_softwares}/multiqc/env



cd ${VC_code}/1.QC_output/

ml miniconda3
conda activate ${multi_qc}
cd ${VC_results}/1.Trim_QC/pre_trim
multiqc ./

cd ${VC_results}/1.Trim_QC/post_trim
multiqc ./

conda deactivate