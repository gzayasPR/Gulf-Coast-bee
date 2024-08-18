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
vcf_file=${PG_data}/variants/${variant_caller}/Hesperapis_oraria.polymorphic.vcf 
meta_data=${PG_data}/variants/Samples.Metadata.csv
name=Hesperapis_oraria
min_dp="4"
min_locus_missing="0.01"
min_indiv_missing="0.01"
min_MAC="2"
prune="TRUE"
mkdir -p ${out_dir}
cd $r_library

ml r/4.4.0
Rscript $PG_code/scripts/Parameters.QC_LD.R ${out_dir} ${vcf_file} ${meta_data} ${name} ${min_dp} ${min_locus_missing} ${min_indiv_missing} ${min_MAC} ${prune}


VCF=${out_dir}/${name}_pruned.vcf.gz
Rscript $PG_code/scripts/PCA.V2.R ${out_dir} ${meta_data} $VCF 