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
source ~/.bashrc

# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $PG_code"
echo "Data Directory: $PG_data"
echo "Results Directory: $PG_results"
# Set the correct environment variables
export PATH=${my_softwares}/bin:$PATH
export LD_LIBRARY_PATH=${my_softwares}/lib:$LD_LIBRARY_PATH
export CFLAGS="-I${my_softwares}/include"
export LDFLAGS="-L${my_softwares}/lib"

variant_caller=GATK
out_dir=${PG_results}/Genetic_Diversity/${variant_caller}
vcf_file=${PG_data}/variants/${variant_caller}/Females_Hesperapis_oraria.vcf 
meta_data=${PG_data}/variants/Samples.Metadata.csv
name=Hesperapis_oraria_Hetero
king=${my_softwares}/king
ngsRelate_dir=${my_softwares}/ngsRelate
mkdir -p ${out_dir}
cd ${out_dir}
rm ${out_dir}/*
ml vcftools
ml bcftools 

min_Q=20
min_meanDP=3 # Setting the minimum read depth to 10 as per the new criteria
max_meanDP=55
max_locus_missing="0.95"
GQ_threshold=20
minDP=10 #Setting the minimum read depth to 10 as per the new criteria
maxDP=100

mkdir -p ${out_dir}
ml vcftools
ml bcftools
bcftools +fixploidy ${vcf_file}  -- -f 2 > ${out_dir}/diploid.vcf

# Separate VCF files by site using metadata
grep "USA-FL_Escambia_co" ${meta_data} | awk -F "," '{print $1}' > ${out_dir}/FL.names
grep "USA-AL_Baldwin_co" ${meta_data} | awk -F "," '{print $1}' > ${out_dir}/AL.names

# Create site-specific VCF files
vcftools --vcf ${out_dir}/diploid.vcf --keep ${out_dir}/FL.names --out ${out_dir}/FL_site --recode --recode-INFO-all
vcftools --vcf ${out_dir}/diploid.vcf --keep ${out_dir}/AL.names --out ${out_dir}/AL_site --recode --recode-INFO-all

# Filter each VCF file
vcftools --vcf ${out_dir}/FL_site.recode.vcf --min-meanDP ${min_meanDP} --max-meanDP ${max_meanDP} \
            --minQ ${min_Q} --max-missing ${max_locus_missing} --out ${out_dir}/FL_site_filtered --recode --recode-INFO-all
vcftools --vcf ${out_dir}/AL_site.recode.vcf --min-meanDP ${min_meanDP} --max-meanDP ${max_meanDP} \
            --minQ ${min_Q} --max-missing ${max_locus_missing} --out ${out_dir}/AL_site_filtered --recode --recode-INFO-all

# Run heterozygosity analysis on filtered VCFs
vcftools --vcf ${out_dir}/FL_site_filtered.recode.vcf --het --out ${out_dir}/FL_site_heterozygosity
vcftools --vcf ${out_dir}/AL_site_filtered.recode.vcf --het --out ${out_dir}/AL_site_heterozygosity

(cat ${out_dir}/AL_site_heterozygosity.het && sed '1d' ${out_dir}/FL_site_heterozygosity.het) > ${out_dir}/${name}.het

bcftools query -l ${out_dir}/FL_site_filtered.recode.vcf > ${out_dir}/FL.pop_info
bcftools query -l ${out_dir}/AL_site_filtered.recode.vcf > ${out_dir}/AL.pop_info


$ngsRelate_dir/ngsRelate/ngsRelate -h ${out_dir}/FL_site_filtered.recode.vcf -n $(wc -l  ${out_dir}/FL.pop_info) -F 1 -T GT -O ${out_dir}/inbreeding.FL
$ngsRelate_dir/ngsRelate/ngsRelate -h ${out_dir}/AL_site_filtered.recode.vcf -n $(wc -l ${out_dir}/AL.pop_info) -F 1 -T GT -O ${out_dir}/inbreeding.AL


( cat ${out_dir}/AL.pop_info && cat ${out_dir}/FL.pop_info ) > ${out_dir}/pop.info
( cat ${out_dir}/inbreeding.AL && sed '1d' ${out_dir}/inbreeding.FL ) > ${out_dir}/inbreeding
cd $r_library 
ml r/4.4.0
Rscript $PG_code/scripts/Hetero.V3.R ${out_dir}/ ${out_dir}/${name}.het ${meta_data}