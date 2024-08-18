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

vcf_dir=${VC_results}/3.5.Intial_calling_V2/
output_dir=${VC_results}/3.5.Intial_calling_V2/Merged_output/
mkdir -p ${output_dir}
# Define the path to the gvcf_list.txt
vcf_list=${output_dir}/vcf_list.txt
#rm -f $vcf_list
# Iterate over each GVCF file and add to the list

#for vcf in $(ls ${vcf_dir}/sorted_genome_intervals_*.vcf); do
#    name=$(basename ${IID} .vcf)
#    echo -e "${vcf}" >> $vcf_list
#done
ml python
ml gatk
gatk  --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=48" GatherVcfs -I $vcf_list -O ${output_dir}/combined.vcf
ml bcftools


bcftools view  ${output_dir}/combined.vcf -H -Oz -o  ${output_dir}/combined.vcf.gz

# Define input and output files
input_vcf="${output_dir}/combined.vcf"
output_vcf="${output_dir}/fixed_input.vcf.gz"


# Validate and preprocess the VCF file using bcftools
bcftools view ${input_vcf} | bcftools norm -m -any -f ${ref_genome} -Oz -o ${output_vcf}

# Index the output VCF file
bcftools index ${output_vcf}
ml vcftools
vcftools --gzvcf ${output_dir}/fixed_input.vcf.gz --min-meanDP 5 --minQ 9  --max-missing 1 --out ${output_dir}/fixed_missing.25 --recode --recode-INFO-all

cd $r_library
name=intial.calling

ml r/4.4.0
Rscript $VC_code/scripts/Filter.Het.R ${output_dir} ${output_dir}/fixed_missing.25.recode.vcf $meta_data $name
