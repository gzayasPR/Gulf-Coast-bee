#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=Variant_%j.out
#SBATCH --error=Variant_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb

source VC_project_env.sh
# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $VC_code"
echo "Data Directory: $VC_data"
echo "Results Directory: $VC_results"

ref_genome=${VC_data}/final_assembly/Hesperapis_oraria_2.curatedScaff.clean_sort_rename_nuc
bam_dir=${VC_results}/3.Cleaning_V2/
known_vcf=${VC_results}/3.5.Intial_calling_V2/Merged_output/intial.calling_no.het.vcf.gz
cd ${VC_data}/final_assembly/
ml python
ml samtools
ml gatk
ml bcftools


# Check if the SAMtools index file exists
#if [ ! -f "${ref_genome}.fasta .fai" ]; then.
#    echo "Index not found for ref genome. Indexing now..."
#    gatk CreateSequenceDictionary -R ${ref_genome}.fasta 
#    samtools faidx ${ref_genome}.fasta 
#else
#    echo "Index already exists. Proceeding with analysis..."
#fi
cd ${VC_code}/4.Calling_output/
ref_genome=${ref_genome}.fasta
output_dir=${VC_results}/4.Calling_V2/
mkdir -p ${output_dir}
mkdir -p ${VC_code}/4.Calling_output/
echo "BLX2737	female	2
BLX2738	female	2
BLX2739	male	1
BLX2740	male	1
BLX2741	male	1
BLX2743	female	2
BLX2744	female	2
BLX2745	male	1
BLX2747	female	2
BLX2748	female	2
BLX2751	male	1
BLX2752	male	1
BLX2753	female	2
BLX2754	female	2
BLX2755	male	1
BLX2756	male	1
BLX2757	male	1
BLX2759	female	2" > ${output_dir}/ploidy.txt

bcftools view ${known_vcf} > ${output_dir}/known.vcf
if [ ! -f "${output_dir}/known.vcf.tbi" ]; then
    gatk IndexFeatureFile -I ${output_dir}/known.vcf
fi
# Iterate over each file ending with "_R1.fastq.gz" in the skimseek directory
for IID in ${bam_dir}/*.marked_duplicates.bam; do
    # Extract the base filename without the suffix "_R1.fastq.gz"
    name=$(basename ${IID} .marked_duplicates.bam )
    echo ${name}
    cd ${VC_code}/4.Calling_output/
    ploidy=$(awk -v name="$name" '$1==name {print $3}' ${output_dir}/ploidy.txt)
    echo $ploidy
    mkdir -p ${output_dir}/GATK
    mkdir -p ${output_dir}/DeepVariant
    #Modify the job-name with the sample name
  #sbatch --job-name="${name}.GATK" --output="${name}.GATK.out" --error="${name}.GATK.err" ../scripts/variant.calling.sh ${name} ${bam_dir} ${ref_genome} ${known_vcf} ${ploidy} ${output_dir}
#  sbatch --job-name="${name}.GATK" --output="${name}.GATK.out" --error="${name}.GATK.err" ../scripts/variant.calling.sh ${name} ${bam_dir} ${ref_genome} ${output_dir}/known.vcf 2 ${output_dir}
    sbatch --job-name="${name}.DeepVariant" --output="${name}.DeepVariant.out" --error="${name}.DeepVariant.err" ../scripts/deep.variant.calling.sh ${name} ${bam_dir} ${ref_genome} ${known_vcf} ${ploidy} ${output_dir}
done

