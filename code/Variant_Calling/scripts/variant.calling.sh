#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --job-name="$1"
#SBATCH --output=$1.out
#SBATCH --error=$1.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=72:00:00  # Adjusted walltime to 24 hours
#SBATCH --ntasks=48  # Using all available cores on a single node
#SBATCH --nodes=1  # Use a single node
#SBATCH --partition=atlas  # Selecting the atlas partition
#SBATCH --qos=normal  # Using the normal QoS


name=$1
bam_dir=$2
ref_genome=$3
known_vcf=$4
ploidy=$5
output_dir=$6/GATK

mkdir -p ${output}
ml gatk
ml samtools
ml picard
ml qualimap

ml miniconda3 
eval "$(conda shell.bash hook)"
input_bam=${bam_dir}/${name}.marked_duplicates.bam

mkdir -p ${output_dir}
# Step 1: Initial Variant Calling
# Step 2: Run BaseRecalibrator
gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=48" BaseRecalibrator \
    -I ${input_bam} \
    -R ${ref_genome} \
    --known-sites ${known_vcf} \
    -O ${output_dir}/recal_data.table.${name}

# Step 3: Apply BQSR
gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=48" ApplyBQSR \
    -R ${ref_genome} \
    -I ${input_bam} \
    -bqsr ${output_dir}/recal_data.table.${name} \
    -O ${output_dir}/${name}.recalibrated.bam

# Step 4: Run HaplotypeCaller
gatk --java-options "-Xmx8g" HaplotypeCaller \
    -R ${ref_genome} \
    -I ${output_dir}/${name}.recalibrated.bam \
    -O ${output_dir}/${name}.g.vcf.gz \
    -A Coverage -A FisherStrand -A StrandOddsRatio \
    -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
    --pcr-indel-model CONSERVATIVE -ploidy ${ploidy} -stand-call-conf 30 \
    -ERC GVCF 