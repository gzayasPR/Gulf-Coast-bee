#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=variants.setup_%j.out
#SBATCH --error=variants.setup_%j.err
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

vcf_file="${PG_data}/variants/Hesperapis_oraria.vcf.gz" 
ml bcftools
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${vcf_file} > ${PG_data}/variants/snp_info.txt
# Extract necessary fields
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${vcf_file} > ${PG_data}/variants/snp_info.txt
# Create new SNP IDs
#awk 'BEGIN{FS=OFS="\t"} {print $1 ":" $2 ":" $3 ":" $4}' ${PG_data}/variants/snp_info.txt > ${PG_data}/variants/snp_ids.txt

# Combine the new SNP IDs with original SNP info
#paste ${PG_data}/variants/snp_ids.txt ${PG_data}/variants/snp_info.txt > ${PG_data}/variants/new_snp_info.txt
# Update the VCF file with new SNP IDs
zcat ${vcf_file} | awk 'BEGIN{FS=OFS="\t"} NR==FNR {a[$1]=$2 ":" $3 ":" $4 ":" $5; next} /^#/ {print; next} {if ($3 == ".") $3=a[$1]; print}' ${PG_data}/variants/new_snp_info.txt > ${PG_data}/variants/output_with_snpid.vcf




