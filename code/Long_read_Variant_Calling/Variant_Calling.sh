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

source LR_VC_project_env.sh
# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $LR_VC_code"
echo "Data Directory: $LR_VC_data"
echo "Results Directory: $LR_VC_results"

bam_dir=${LR_VC_results}/Alignment/
ml python
ml samtools
ml gatk
ml bcftools
output_dir=${LR_VC_results}/Variant.Calling/
mkdir -p ${output_dir}
mkdir -p ${LR_VC_code}/Calling_output/
cd ${LR_VC_code}/Calling_output/
# Process the short reads metadata
# Iterate over each file ending with "_R1.fastq.gz" in the skimseek directory
for IID in ${bam_dir}/*.bam; do
    # Extract the base filename without the suffix "_R1.fastq.gz"
    name=$(basename ${IID} .bam )
    echo ${name}
    cd  ${LR_VC_code}/Calling_output/
    mkdir -p ${output_dir}/DeepVariant
    #Modify the job-name with the sample name
    #sbatch --job-name="${name}.DeepVariant" --output="${name}.DeepVariant.out" --error="${name}.DeepVariant.err" ../scripts/deep.variant.calling.sh ${name} ${bam_dir} ${ref_genome} ${known_LR_VCf} ${ploidy} ${output_dir} ${my_softwares}
    sbatch --job-name="${name}.DeepVariant" --output="${name}.DeepVariant.out" --error="${name}.DeepVariant.err" ../scripts/deep.variant.calling.sh ${name} ${bam_dir} ${ref_genome} 2 ${output_dir} ${my_softwares}
done

