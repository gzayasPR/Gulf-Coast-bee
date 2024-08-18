#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=Final_calling_%j.out
#SBATCH --error=Final_calling_%j.err
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

clean_dir=${VC_results}/4.Calling/
output_dir=${VC_results}/4.5.V2.Final_Calling/
mkdir -p ${output_dir}
# Define the path to the gvcf_list.txt
gvcf_list=${output_dir}/gvcf_list.txt

# Remove existing gvcf_list.txt if it exists
#rm -f $gvcf_list
ml bcftools
# Iterate over each GVCF file and add to the list
#for IID in ${clean_dir}/*.g.vcf.gz; do
#  name=$(basename ${IID} .g.vcf.gz)
#   echo -e "${name}\t${IID}" >> $gvcf_list
#   bcftools index -f ${IID}
#done

# Generate the interval file (genome_intervals.bed or genome_intervals.list)
#awk '{print $1"\t0\t"$2}' ${ref_genome}.fai > ${output_dir}/genome_intervals.bed


ml gatk
#rm -R ${output_dir}/genomicsdb 
# Combine GVCF files using GenomicsDBImport
#gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=48" GenomicsDBImport \
#    --genomicsdb-workspace-path ${output_dir}/genomicsdb \
#    --batch-size 50 \
#    --L ${output_dir}/genome_intervals.bed \
#    --sample-name-map $gvcf_list \
#    --merge-input-intervals \
#    --merge-contigs-into-num-partitions 1

# Joint Genotyping
gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=48" GenotypeGVCFs \
    -R ${ref_genome} \
    -V gendb://${output_dir}/genomicsdb \
    --include-non-variant-sites \
    -O ${output_dir}/combined_variants.vcf
   
  

#gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=48" VariantFiltration \
#   -R ${ref_genome} \
#    -V ${output_dir}/combined_variants.vcf \
#    -O ${output_dir}/filtered_variants.vcf \
#     --filter-name "QD2" --filter-expression "QD < 1.0" \
#    --filter-name "FS60" --filter-expression "FS > 100.0" \
#    --filter-name "MQ40" --filter-expression "MQ < 30.0" \
#    --filter-name "MQRankSum-12.5" --filter-expression "MQRankSum < -15.0" \
#    --filter-name "ReadPosRankSum-8" --filter-expression "ReadPosRankSum <  -10.0 "

# Index the VCF files
ml bcftools
#bcftools view  ${output_dir}/filtered_variants.vcf -Oz -o  ${output_dir}/filtered_variants.vcf.gz
#bcftools index -f ${output_dir}/filtered_variants.vcf.gz

# Select only passing variants
#gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=48" SelectVariants \
#    -R ${ref_genome} \
#    -V ${output_dir}/filtered_variants.vcf  \
#    -O ${output_dir}/high_quality_variants.vcf \
#    --set-filtered-gt-to-nocall \
#    --exclude-filtered



#bcftools view  ${output_dir}/high_quality_variants.vcf -Oz -o  ${output_dir}/high_quality_variants.vcf.gz
#bcftools index -f ${output_dir}/high_quality_variants.vcf.gz
#bcftools +setGT ${output_dir}/high_quality_variants.vcf.gz -- -t q -n . -i 'FMT/DP=0' >  ${output_dir}/high_quality_variants.na.vcf
#bcftools view   ${output_dir}/high_quality_variants.na.vcf -Oz -o   ${output_dir}/high_quality_variants.na.vcf.gz
#bcftools index -f ${output_dir}/high_quality_variants.na.vcf.gz
#bcftools view -H ${output_dir}/high_quality_variants.na.vcf.gz | head -20


#ml plink2
# Define input and output files
#input_vcf="${output_dir}/high_quality_variants.na.vcf.gz"
#output_vcf="${output_dir}/fixed_input.vcf.gz"


# Validate and preprocess the VCF file using bcftools
#bcftools view ${input_vcf} | bcftools norm -m -any -f ${ref_genome} -Oz -o ${output_vcf}

# Index the output VCF file
#bcftools index ${output_vcf}

# Check the number of variants in the fixed VCF file
#echo "Number of variants in the fixed VCF:"
#bcftools view -H ${output_vcf} | wc -l
#ml plink2
#plink2 --allow-extra-chr --vcf ${output_dir}/fixed_input.vcf.gz --freq --out allele_freq
#awk '$5 == 0 || $5 == 1' ${output_dir}/allele_freq.afreq > ${output_dir}/monomorphic_sites.txt

# Set the total number of samples in your VCF file
#total_samples=18  # Replace with the actual number of samples

#for missing in 9 7 4 2 1 0; do
    # Convert the missing count to a proportion for --geno
#    geno_threshold=$(echo "scale=10; ${missing} / ${total_samples}" | bc)

#    echo "Filtering with --geno threshold: ${geno_threshold} which is $missing"
    
    # Run plink2 with the calculated geno threshold
#    plink2 --allow-extra-chr --vcf ${output_dir}/fixed_input.vcf.gz \
#           --geno ${geno_threshold} \
#           --recode vcf \
#           --out ${output_dir}/plink.filtered_variants_${missing}
#done

#plink2 --allow-extra-chr --vcf ${output_dir}/plink.filtered_variants_0.vcf --ref-from-fa ${ref_genome} -make-bed -out ${output_dir}/allele_freq_0
#plink2 --allow-extra-chr --bfile ${output_dir}/allele_freq_0 -make-bed --freq --out ${output_dir}/allele_freq_0
#awk '$5 == 0 || $5 == 1' ${output_dir}/allele_freq.afreq_0 > ${output_dir}/monomorphic_sites_0.txt
