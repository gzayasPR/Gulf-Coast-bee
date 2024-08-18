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
echo "Results Directory: $PG_results"

variant_caller=GATK
out_dir=${PG_results}/Genetic_Diversity/${variant_caller}
vcf_file=${PG_data}/variants/${variant_caller}/Hesperapis_oraria.vcf 
meta_data=${PG_data}/variants/Samples.Metadata.csv
name=Hesperapis_oraria_Hetero
king=${my_softwares}/king

mkdir -p ${out_dir}
cd ${out_dir}
rm ${out_dir}/*
ml vcftools
ml bcftools 

min_Q=20
min_meanDP=3
max_meanDP=55
max_locus_missing="0.95"
mkdir -p ${out_dir}
ml vcftools
vcftools --vcf ${vcf_file} --max-missing ${max_locus_missing} --min-meanDP ${min_meanDP} --max-meanDP ${max_meanDP} --minQ ${min_Q} --out ${out_dir}/fixed_missing --recode --recode-INFO-all
vcftools --vcf ${out_dir}/fixed_missing.recode.vcf --het --out ${out_dir}/${name}
cd $r_library 
ml r/4.4.0
Rscript $PG_code/scripts/Hetero.V3.R ${out_dir}/ ${out_dir}/${name}.het ${meta_data}



ml plink2
plink2 --vcf  ${vcf_file} --allow-extra-chr  --vcf-half-call h  --make-king triangle --make-bed --out ${out_dir}/filtered_data
awk 'BEGIN{OFS="\t"} { $1=NR; print }' ${out_dir}/filtered_data.bim > ${out_dir}/filtered_data_autosomal.bim

cp filtered_data.bed filtered_data_autosomal.bed
cp filtered_data.fam filtered_data_autosomal.fam
$king -b ${out_dir}/filtered_data_autosomal.bed --kinship --sexchr $(tail -1 ${out_dir}/filtered_data_autosomal.bim  | cut -f1 )  --prefix ${out_dir}/kinship_coeff
cat kinship_coeff.kin 

