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

bam_dir=${VC_results}/3.Cleaning/
known_vcf=${VC_results}/3.5.Intial_calling/Merged_output/fixed_missing.25.recode.vcf
cd ${VC_data}/final_assembly/
ml python
ml samtools
ml gatk
ml bcftools

cd ${VC_code}/4.Calling_output/
output_dir=${VC_results}/4.Calling/
mkdir -p ${output_dir}
mkdir -p ${VC_code}/4.Calling_output/
# Process the short reads metadata
awk -F',' '
NR==1 {for (i=1; i<=NF; i++) {if ($i=="BioSample") b=i; if ($i=="Sex") f1=i}}
NR>1 {
    ploidy = ($f1 == "female") ? 2 : 1
    print $b","$f1","ploidy
}' "$meta_data" > ${output_dir}/ploidy.txt

bcftools view ${known_vcf} > ${output_dir}/known.vcf
if [ ! -f "${output_dir}/known.vcf.tbi" ]; then
    gatk IndexFeatureFile -I ${output_dir}/known.vcf
fi
known_vcf=${output_dir}/known.vcf
# Iterate over each file ending with "_R1.fastq.gz" in the skimseek directory
for IID in ${bam_dir}/*.marked_duplicates.bam; do
    # Extract the base filename without the suffix "_R1.fastq.gz"
    name=$(basename ${IID} .marked_duplicates.bam )
    echo ${name}
    cd ${VC_code}/4.Calling_output/
    ploidy=$(awk -F',' -v name="$name" '$1==name {print $3}' ${output_dir}/ploidy.txt)
    if [ -z "$ploidy" ]; then
        echo "Ploidy not found for ${name}. Skipping this sample."
        continue
    fi
    echo "Ploidy for ${name}: $ploidy"
    mkdir -p ${output_dir}/GATK
    mkdir -p ${output_dir}/DeepVariant
    #Modify the job-name with the sample name
  sbatch --job-name="${name}.GATK" --output="${name}.GATK.out" --error="${name}.GATK.err" ../scripts/variant.calling.sh ${name} ${bam_dir} ${ref_genome} ${known_vcf} ${ploidy} ${output_dir} ${my_softwares}
    #sbatch --job-name="${name}.DeepVariant" --output="${name}.DeepVariant.out" --error="${name}.DeepVariant.err" ../scripts/deep.variant.calling.sh ${name} ${bam_dir} ${ref_genome} ${known_vcf} ${ploidy} ${output_dir} ${my_softwares}
done

