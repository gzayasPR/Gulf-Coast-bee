#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=MergingVCF_%j.out
#SBATCH --error=MergingVCF_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=168:00:00
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
vcf_dir=${VC_results}/4.5.Final_Calling/
output_dir=${VC_results}/5.Merge_Calling/GATK/
mkdir -p ${output_dir}

# Define the path to the vcf_list.txt
vcf_list=${output_dir}/vcf_list.txt
rm -f $vcf_list

# Iterate over each VCF file and add to the list
for vcf in ${vcf_dir}/sorted_genome_intervals_*.vcf; do
    echo -e "${vcf}" >> $vcf_list
done

module load gatk
module load python
rm -f ${output_dir}/combined.vcf
gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=48" GatherVcfs -I $vcf_list -O ${output_dir}/combined.vcf

module load bcftools
bcftools view ${output_dir}/combined.vcf -H -Oz -o ${output_dir}/combined.vcf.gz

# Define input and output files
input_vcf="${output_dir}/combined.vcf.gz"
output_vcf="${output_dir}/fixed_input.vcf.gz"

# Validate and preprocess the VCF file using bcftools
bcftools view ${input_vcf} | bcftools norm -m -any -f ${ref_genome} -Oz -o ${output_vcf}

# Index the output VCF file
bcftools index ${output_vcf}

module load vcftools
vcftools --gzvcf ${output_vcf} --max-missing 0.01 --out ${output_dir}/fixed_missing --recode --recode-INFO-all

mkdir -p ${VC_data}/variants/GATK
cp ${output_dir}/fixed_missing.recode.vcf ${VC_data}/variants/GATK/temp.vcf

module load bcftools
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' ${VC_data}/variants/GATK/temp.vcf > ${VC_data}/variants/GATK/Hesperapis_oraria.vcf

# Filter for minor allele count of 1 or greater
bcftools view -i 'MAC>=1' ${VC_data}/variants/GATK/Hesperapis_oraria.vcf -Oz -o ${VC_data}/variants/GATK/Hesperapis_oraria.polymorphic.vcf

# Generate females-specific VCF
grep "female" $meta_data | awk -F "," '{print $1}' > ${VC_data}/variants/GATK/Females.ID
bcftools view -S ${VC_data}/variants/GATK/Females.ID ${VC_data}/variants/GATK/Hesperapis_oraria.polymorphic.vcf -o ${VC_data}/variants/GATK/Females_Hesperapis_oraria.vcf

