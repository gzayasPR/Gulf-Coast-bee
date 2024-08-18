#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=Taj_D_%j.out
#SBATCH --error=Taj_D_%j.err
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
out_dir=${PG_results}/Taj_D/DeepVariant/
vcf_file=${PG_data}/variants/DeepVariant/Hesperapis_oraria_no.het.vcf.gz
meta_data=${PG_data}/variants/Samples.Metadata.csv
name=Hesperapis_oraria
min_dp="4"
min_locus_missing="0.15"
min_indiv_missing="0.15"
min_MAC="0"
prune="FALSE"
mkdir -p ${out_dir}
cd $r_library

ml r/4.4.0
Rscript $PG_code/scripts/Parameters.QC_LD.R ${out_dir} ${vcf_file} ${meta_data} ${name} ${min_dp} ${min_locus_missing} ${min_indiv_missing} ${min_MAC} ${prune}

cd ${out_dir}
ml vcftools
vcftools --gzvcf ${out_dir}/${name}_filtered.vcf.gz --TajimaD 100000 --out ${out_dir}/${name}
vcftools  --gzvcf ${out_dir}/${name}_filtered.vcf.gz --window-pi 100000 --out ${out_dir}/${name}
cd ${PG_code}

cd $r_library

Rscript $PG_code/scripts/TajimasD.R ${out_dir} ${meta_data}  ${out_dir}/${name}.Tajima.D
#Rscript $PG_code/scripts/TajimasD.R ${out_dir} ${meta_data}  ${out_dir}/${name}.Tajima.D
cd ${PG_code}



out_dir=${PG_results}/Taj_D/GATK/
vcf_file=${PG_data}/variants/Hesperapis_oraria.vcf
meta_data=${PG_data}/variants/Samples.Metadata.csv
name=Hesperapis_oraria
min_dp="4"
min_locus_missing="0.15"
min_indiv_missing="0.15"
min_MAC="0"
prune="FALSE"
mkdir -p ${out_dir}
cd $r_library

ml r/4.4.0
Rscript $PG_code/scripts/Parameters.QC_LD.R ${out_dir} ${vcf_file} ${meta_data} ${name} ${min_dp} ${min_locus_missing} ${min_indiv_missing} ${min_MAC} ${prune}

cd ${out_dir}
ml bcftools
bcftools +fixploidy ${out_dir}/${name}_filtered.vcf.gz -- -f 2 > ${out_dir}/${name}_diploid.vcf


ml vcftools
vcftools --vcf ${out_dir}/${name}_diploid.vcf --TajimaD 100000 --out ${out_dir}/${name}
vcftools  --vcf ${out_dir}/${name}_diploid.vcf --window-pi 100000 --out ${out_dir}/${name}
cd ${PG_code}

cd $r_library

Rscript $PG_code/scripts/TajimasD.R ${out_dir} ${meta_data}  ${out_dir}/${name}.Tajima.D
#Rscript $PG_code/scripts/TajimasD.R ${out_dir} ${meta_data}  ${out_dir}/${name}.Tajima.D
cd ${PG_code}