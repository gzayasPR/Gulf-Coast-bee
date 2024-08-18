#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --job-name="runClean_$1"
#SBATCH --output=runClean_$1_%j.out
#SBATCH --error=runClean_$1_%j.err
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
ploidy=$4
output_dir=$5
my_softwares=$6

ml python
ml gatk
ml samtools
ml picard
echo "starting rg"
picard AddOrReplaceReadGroups \
    I= ${bam_dir}/${name}.sorted_aligned_reads.bam  \
    O=${output_dir}/${name}.rg_added.bam \
    RGID=1 \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=${name}
echo "finished rg"    
echo "starting Marking Duplicates"
gatk --java-options "-Xmx4g -XX:+UseParallelGC -XX:ParallelGCThreads=48" MarkDuplicates -I ${output_dir}/${name}.rg_added.bam -O ${output_dir}/${name}.marked_duplicates.bam -M ${output_dir}/${name}.marked_dup_metrics.txt
samtools index  -@ 48 ${output_dir}/${name}.marked_duplicates.bam
echo "Finished Marking Duplicates"
ml miniconda3 
eval "$(conda shell.bash hook)"
conda activate ${my_softwares}/qualimap/env
export JAVA_OPTS="-Xms4G -Xmx32G" 
echo "starting qualimap"
qualimap bamqc -bam ${output_dir}/${name}.marked_duplicates.bam 
echo "finished qualimap"
# Step 1: Initial Variant Calling
echo "starting haplotypecaller"
gatk  --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=48" HaplotypeCaller \
  -R ${ref_genome} \
  -I ${output_dir}/${name}.marked_duplicates.bam  \
  -O ${output_dir}/${name}.g.vcf.gz \
  -A Coverage -A FisherStrand -A StrandOddsRatio \
  -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
  --pcr-indel-model CONSERVATIVE -ploidy ${ploidy} -stand-call-conf 30 \
  -ERC GVCF 
echo "Finished haplotypecaller"