#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=GD_%j.out
#SBATCH --error=GD_%j.err
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

vairant_caller=DeepVariant
out_dir=${PG_results}/Genetic_Diversity/${vairant_caller}
vcf_file=${PG_data}/variants/${vairant_caller}/Females_Hesperapis_oraria.vcf
meta_data=${PG_data}/variants/Samples.Metadata.csv
name=Hesperapis_oraria_Hetero
min_dp="5"
min_locus_missing="0.25"
min_indiv_missing="0.25"
min_MAC="0"
prune="FALSE"
mkdir -p ${out_dir}
cd ${out_dir}
ml vcftools
ml bcftools 
bcftools +setGT --threads 48 ${vcf_file} -- -t q -n . -i "FMT/DP<$min_dp" > ${out_dir}/temp.vcf
ml vcftools
vcftools --vcf ${out_dir}/temp.vcf --max-missing $min_locus_missing --out ${out_dir}/temp2.vcf --recode --recode-INFO-all

cd $r_library

ml r/4.4.0
Rscript $PG_code/scripts/Parameters.QC_LD.Females.R ${out_dir} ${out_dir}/temp2.vcf.recode.vcf ${meta_data} ${name} ${min_dp} ${min_locus_missing} ${min_indiv_missing} ${min_MAC} ${prune}
cd ${out_dir}
mkdir -p ${out_dir}
cd $r_library

ml r/4.4.0
Rscript $PG_code/scripts/Hetero.V2.R ${out_dir}/ ${out_dir}/temp2.vcf.recode.vcf ${meta_data}


#out_dir=${PG_results}/PCA/DeepVariant
#mkdir -p ${out_dir}
#cd $r_library

#ml r/4.4.0
#Rscript $PG_code/scripts/DV.QC_LD.R
#Rscript $PG_code/scripts/DV.PCA.R