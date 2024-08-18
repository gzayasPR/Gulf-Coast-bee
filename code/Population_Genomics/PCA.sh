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
variant_caller=GATK
out_dir=${PG_results}/PCA/${variant_caller}/
vcf_file=${PG_data}/variants/${variant_caller}/Hesperapis_oraria.vcf 
meta_data=${PG_data}/variants/Samples.Metadata.csv
name=Hesperapis_oraria
min_Q=20
min_meanDP=3
max_meanDP=55
max_locus_missing="0.95"
MAC=3
mkdir -p ${out_dir}
ml vcftools
ml plink2
vcftools --vcf ${vcf_file} --max-missing ${max_locus_missing} --mac ${MAC} --min-meanDP ${min_meanDP} --max-meanDP ${max_meanDP} --minQ ${min_Q} --out ${out_dir}/fixed_missing --recode --recode-INFO-all
plink2 --vcf ${out_dir}/fixed_missing.recode.vcf --vcf-half-call h  --allow-extra-chr -recode vcf -out ${out_dir}/fixed

cd $r_library
VCF=${out_dir}/fixed.vcf


cd $r_library
ml r/4.4.0
Rscript $PG_code/scripts/PCA.V2.R ${out_dir} ${meta_data} $VCF 




cd $PG_code

variant_caller=Deep_Variant
out_dir=${PG_results}/PCA/${variant_caller}/
vcf_file=${PG_data}/variants/${variant_caller}/Hesperapis_oraria.vcf 
meta_data=${PG_data}/variants/Samples.Metadata.csv
name=Hesperapis_oraria
min_Q=20
min_meanDP=3
max_meanDP=55
max_locus_missing="0.95"
MAC=3
mkdir -p ${out_dir}
ml vcftools
ml plink2
vcftools --vcf ${vcf_file} --max-missing ${max_locus_missing} --mac ${MAC} --min-meanDP ${min_meanDP} --max-meanDP ${max_meanDP} --minQ ${min_Q} --out ${out_dir}/fixed_missing --recode --recode-INFO-all
plink2 --vcf ${out_dir}/fixed_missing.recode.vcf --vcf-half-call h  --allow-extra-chr -recode vcf -out ${out_dir}/fixed

cd $r_library
VCF=${out_dir}/fixed.vcf


cd $r_library
ml r/4.4.0
Rscript $PG_code/scripts/PCA.V2.R ${out_dir} ${meta_data} $VCF 
