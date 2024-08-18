#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=PCA_%j.out
#SBATCH --error=PCA_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb

source PG_project_env.sh
# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $PG_code"
echo "Data Directory: $PG_data"
echo "Results Directory: $PG_results"Í

mkdir -p ${PG_results}/PCA/
out_dir=${PG_results}/PCA/Compare_PCA
mkdir -p ${out_dir}
PCA_Scores=${PG_code}/PCA_Scores.csv

cd $r_library

ml r/4.4.0
Rscript ${PG_code}/scripts/compare.PCA_scores.R ${out_dir} ${PCA_Scores}


