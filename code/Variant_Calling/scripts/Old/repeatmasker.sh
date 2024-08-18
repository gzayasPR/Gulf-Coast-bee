#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=runClean_%j.out
#SBATCH --error=runClean_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb

source VC_project_env.sh
# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $VC_code"
echo "Data Directory: $VC_data"
echo "Results Directory: $VC_results"

ml miniconda3
eval "$(conda shell.bash hook)"

conda activate ${VC_conda}
#conda install -c bioconda repeatmasker
cd  $VC_data/final_assembly
RepeatMasker -species bees Hesperapis_oraria_2.curatedScaff.clean_sort_rename_nuc.fasta
