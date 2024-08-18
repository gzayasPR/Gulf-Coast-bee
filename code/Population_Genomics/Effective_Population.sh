#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=Ne_%j.out
#SBATCH --error=Ne_%j.err
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
echo "Results Directory: $PG_results"

variant_caller=DeepVariant
out_dir=${PG_results}/Genetic_Diversity/${variant_caller}
vcf_file=${PG_data}/variants/${variant_caller}/Hesperapis_oraria.vcf 
meta_data=${PG_data}/variants/Samples.Metadata.csv
LDne=${my_softwares}/ldne/
plink_19=/project/beenome100/conda_envs/plink1.9/
mkdir -p ${out_dir}
cd ${out_dir}
rm ${out_dir}/*
ml vcftools
ml bcftools 
ml plink2/2.00a4.3 

min_Q=20
min_meanDP=3
max_meanDP=55
max_locus_missing="0.95"
mkdir -p ${out_dir}
vcftools --vcf ${vcf_file} --max-missing ${max_locus_missing} --min-meanDP ${min_meanDP} --max-meanDP ${max_meanDP} --minQ ${min_Q} --out ${out_dir}/fixed_missing --recode --recode-INFO-all
ml miniconda3
conda activate ${plink_19}
plink --vcf ${out_dir}/fixed_missing.recode.vcf  -allow-extra-chr --vcf-half-call h --r2 --ld-window-r2 0.2 --out ${out_dir}/ld_output
chmod +x ${LDne}/LDNe


${LDne}/LDNe ${out_dir}/ld_output.ld --output ${out_dir}/ldne_results


