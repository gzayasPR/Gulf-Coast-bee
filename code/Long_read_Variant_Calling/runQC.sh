#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=runQC_%j.out
#SBATCH --error=runQC_%j.err
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
reads_dir=/project/90daydata/beenome100/hesperapis_oraria_genomics/pacbio_filtered_reads/


cd ${LR_VC_results}/
mkdir -p ${LR_VC_results}/Quality_Control/
output_dir=${LR_VC_results}/Quality_Control/
mkdir -p ${LR_VC_code}/Quality_Control_output/
# Iterate over each file ending with "_R1.fastq.gz" in trimmed dir
for IID in ${reads_dir}*.hifi_reads.fcsfilt.fastq.gz; do
    # Extract the base filename without the suffix "_R1.fastq.gz"
    name=$(basename ${IID} .hifi_reads.fcsfilt.fastq.gz)
    echo ${name}
    cd ${LR_VC_code}/Quality_Control_output/
    #Modify the job-name with the sample name
    reads=${IID}
  sbatch --job-name="runAlign_${name}" --output="runAlign_${name}.out" --error="runAlign_${name}.err" ../scripts/quality_control.sh ${name} ${reads}  ${ref_genome} ${output_dir}
done

