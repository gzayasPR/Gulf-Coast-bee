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
bam_dir=${VC_results}/3.Cleaning/
known_vcf=${VC_results}/3.5.Intial_calling/Merged_output/fixed_input.vcf.gz
cd ${VC_data}/final_assembly/
ml samtools
ml gatk
# Check if the SAMtools index file exists
#if [ ! -f "${ref_genome}.fasta .fai" ]; then.
#    echo "Index not found for ref genome. Indexing now..."
#    gatk CreateSequenceDictionary -R ${ref_genome}.fasta 
#    samtools faidx ${ref_genome}.fasta 
#else
#    echo "Index already exists. Proceeding with analysis..."
#fi
if [ ! -f "${known_vcf}.tbi" ]; then
    gatk IndexFeatureFile -I ${known_vcf}
fi
cd ${VC_code}/4.Calling_output/
ref_genome=${ref_genome}.fasta
mkdir -p ${VC_results}/4.Calling
output_dir=${VC_results}/4.Calling
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


# Iterate over each file ending with "_R1.fastq.gz" in the skimseek directory
for IID in ${bam_dir}/*.marked_duplicates.bam; do
    # Extract the base filename without the suffix "_R1.fastq.gz"
    name=$(basename ${IID} .marked_duplicates.bam )
    echo ${name}
    cd ${VC_code}/4.Calling_output//
    ploidy=$(awk -v name="$name" '$1==name {print $3}' ${output_dir}/ploidy.txt)
    echo $ploidy
    #Modify the job-name with the sample name
  #sbatch --job-name="${name}" --output="${name}.out" --error="${name}.err" ../scripts/variant.calling.sh ${name} ${bam_dir} ${ref_genome} ${known_vcf} ${ploidy} ${output_dir}
done

