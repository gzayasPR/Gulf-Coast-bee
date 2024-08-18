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
out_dir=${PG_results}/ROH/DeepVariant/
mkdir -p ${out_dir}
cd $r_library

ml r/4.4.0
Rscript $PG_code/scripts/Females_only.QC.LD.R


vcf_file=${PG_data}/variants/DeepVariant/Females_Hesperapis_oraria_ROH.vcf
eval "$(conda shell.bash hook)"
conda activate /project/beenome100/conda_envs/plink1.9
cd $out_dir

plink --allow-extra-chr -vcf ${vcf_file} --vcf-half-call h --recode --out ${out_dir}/females
plink --allow-extra-chr -vcf ${vcf_file} --vcf-half-call h  --homozyg --homozyg-window-snp 20 --homozyg-snp 15 --homozyg-window-missing 3 --homozyg-kb 10 --homozyg-density 1000 
cp ${out_dir}/females.ped ${out_dir}/females2.ped
cd $r_library


Rscript $PG_code/scripts/ROH.R

#cut -f 1 ${out_dir}/females.map | sed -e 's/[^0-9]//g' > ${out_dir}/map2

#cut -f 2- ${out_dir}/females.map > ${out_dir}/map3
#paste ${out_dir}/map2 ${out_dir}/map3 > ${out_dir}/females2.map
#cd $r_library
#Rscript $PG_code/scripts/ROHs.GATK.R
