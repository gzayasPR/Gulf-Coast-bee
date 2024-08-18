#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=runClean_%j.out
#SBATCH --error=runClean_%j.err
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

bam_dir=${VC_results}/2.Alignment/


ml samtools
ml gatk
# Check if the SAMtools index file exists
if [ ! -f "${ref_genome}.fai" ]; then
    echo "Index not found for ref genome. Indexing now..."
    gatk CreateSequenceDictionary -R ${ref_genome}
    samtools faidx ${ref_genome}
else
    echo "Index already exists. Proceeding with analysis..."
fi

cd ${VC_results}/
output_dir=${VC_results}/3.Cleaning/
mkdir -p ${output_dir}
mkdir -p ${VC_code}/3.Cleaning_output/

# Process the short reads metadata
awk -F',' '
NR==1 {for (i=1; i<=NF; i++) {if ($i=="BioSample") b=i; if ($i=="Sex") f1=i}}
NR>1 {
    ploidy = ($f1 == "female") ? 2 : 1
    print $b","$f1","ploidy
}' "$meta_data" > ${output_dir}/ploidy.txt

# Iterate over each file ending with "_R1.fastq.gz" in the skimseek directory
for IID in ${bam_dir}*.sorted_aligned_reads.bam; do
    # Extract the base filename without the suffix ".sorted_aligned_reads.bam"
    name=$(basename ${IID} .sorted_aligned_reads.bam)
    echo "Processing sample: ${name}"
    cd ${VC_code}/3.Cleaning_output/
    # Retrieve ploidy for the current sample
    ploidy=$(awk -F',' -v name="$name" '$1==name {print $3}' ${output_dir}/ploidy.txt)
    if [ -z "$ploidy" ]; then
        echo "Ploidy not found for ${name}. Skipping this sample."
        continue
    fi
    echo "Ploidy for ${name}: $ploidy"
    
    # Modify the job-name with the sample name and submit the job
     sbatch --job-name="runClean_${name}" --output="runClean_${name}.out" --error="runClean_${name}.err" ../scripts/cleaning.bam.sh ${name} ${bam_dir} ${ref_genome} ${ploidy} ${output_dir} ${my_softwares}
done

