#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=3.5.Intial_calling_%j.out
#SBATCH --error=3.5.Intial_calling_%j.err
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

ref_genome=${VC_data}/final_assembly/Hesperapis_oraria_2.curatedScaff.clean_sort_rename_nuc.fasta

output_dir=${VC_results}/GC_Content
clean_dir=${VC_results}/3.Cleaning/
mkdir -p ${output_dir}

cd ${output_dir}
ml bedtools2 
# Using BEDTools to get GC content
bedtools makewindows -g ${ref_genome}.fai -w 100 > ${output_dir}/windows.bed
bedtools nuc -fi ${ref_genome} -bed ${output_dir}/windows.bed > ${output_dir}/gc_content.txt


bedtools coverage -a ${output_dir}/windows.bed -b ${clean_dir} > coverage.txt