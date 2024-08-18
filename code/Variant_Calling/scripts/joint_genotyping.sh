#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --job-name="chunk_$1"
#SBATCH --output=chunk_$1_%j.out
#SBATCH --error=chunk_$1_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=120:00:00  # Adjusted walltime to 24 hours
#SBATCH --ntasks=48  # Using all available cores on a single node
#SBATCH --nodes=1  # Use a single node
#SBATCH --partition=atlas  # Selecting the atlas partition
#SBATCH --qos=normal  # Using the normal QoS

name=$1
chunk_output_dir=$2
chunk=$3 
gvcf_list=$4 
ref_genome=$5
output_dir=$6
paths=$7

echo "Parameters used "
echo "Chunk = $name 
Chunk Directory = $chunk_output_dir
Reference Genome = $ref_genome
Output Directory = $output_dir"
echo "Bed Intervals"
cat ${chunk}
echo "Samples used"
cat $gvcf_list
source $paths
# Load necessary modules
ml python
ml gatk
echo "Starting GenomicsDBImport"
rm -rf $chunk_output_dir
gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=8" GenomicsDBImport \
            --genomicsdb-workspace-path $chunk_output_dir \
            -R ${ref_genome} \
            --batch-size 50 \
            --L $chunk \
            --sample-name-map $gvcf_list 
echo "Finished GenomicsDBImport"
# Joint Genotyping
echo "Starting GenotypeGVCFs"
rm ${output_dir}/${name}.vcf
gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=8" GenotypeGVCFs \
       -R ${ref_genome} \
       -V gendb://$chunk_output_dir \
        --L $chunk \
        --include-non-variant-sites \
       -O ${output_dir}/${name}.vcf
echo "Finished GenotypeGVCFs"
echo "Starting VariantFiltration"
rm ${output_dir}/filtered_${name}.vcf 
gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=8" VariantFiltration \
   -R ${ref_genome} \
    -V ${output_dir}/${name}.vcf \
    -O ${output_dir}/filtered_${name}.vcf \
#   --filter-expression "QUAL < 30" --filter-name "QUAL_FAIL" --filter-expression "DP < 11" --filter-name "DP_FAIL" \
   --filter-expression "QD < 2.0 || MQ < 40.0 || MQRankSum < -12.4 || ReadPosRankSum < -8.0 || FS > 60.0 || SOR > 3.0" --filter-name  "SNP_HARD_FAIL"
echo "Finished VariantFiltration"
ml bcftools
bcftools view --threads 48 ${output_dir}/filtered_${name}.vcf -Oz -o ${output_dir}/filtered_${name}.vcf.gz
bcftools index --threads 48 -f ${output_dir}/filtered_${name}.vcf.gz

# Select only passing variants
echo "Starting SelectVariant"
gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=8" SelectVariants \
    -R ${ref_genome} \
    -V ${output_dir}/filtered_${name}.vcf  \
    -O ${output_dir}/high_quality_${name}.vcf  \
    --exclude-filtered
echo "Finished SelectVariant"
bcftools +setGT --threads 48 ${output_dir}/high_quality_${name}.vcf -- -t q -n . -i 'FMT/DP=0' > ${output_dir}/high_quality_variants.na.${name}.vcf
bcftools view --threads 48 -i 'F_MISSING<0.99' ${output_dir}/high_quality_variants.na.${name}.vcf -o ${output_dir}/high_quality_variants.no.na.${name}.vcf

echo "Starting SortVcf"
gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=8" SortVcf \
    -I ${output_dir}/high_quality_variants.no.na.${name}.vcf \
    -O ${output_dir}/sorted_${name}.vcf
echo "Finished SortVcf"
echo "$name processed final vcf is at ${output_dir}/sorted_${name}.vcf"






