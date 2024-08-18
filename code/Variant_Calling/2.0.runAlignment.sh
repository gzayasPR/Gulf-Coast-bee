#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=runAlign_%j.out
#SBATCH --error=runAlign_%j.err
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
trimmed_dir=${VC_results}/1.Trim_QC/post_trim/

ml bwa
ml samtools
ml picard
ml gatk

# Check if the BWA index files exist
if [ ! -f "${ref_genome}.bwt" ]; then
    echo "Start Indexing Genome with BWA"
    bwa index -a bwtsw ${ref_genome}
    echo "Finished Indexing Genome with BWA"
else
    echo "BWA index files already exist, skipping indexing."
fi

# Check if the samtools index exists
if [ ! -f "${ref_genome}.fai" ]; then
    echo "Creating FAI index with samtools"
    samtools faidx ${ref_genome}
else
    echo "FAI index already exists, skipping creation."
fi

# Check if the Picard sequence dictionary exists
if [ ! -f "${ref_genome%.fasta}.dict" ]; then
    echo "Creating sequence dictionary with Picard"
    picard CreateSequenceDictionary -R ${ref_genome}
else
    echo "Sequence dictionary already exists, skipping creation."
fi


cd ${VC_results}/
mkdir -p ${VC_results}/2.Alignment
output_dir=${VC_results}/2.Alignment
mkdir -p ${VC_code}/2.Alignment_output/
# Iterate over each file ending with "_R1.fastq.gz" in trimmed dir
for IID in ${trimmed_dir}*_val_1.fq; do
    # Extract the base filename without the suffix "_R1.fastq.gz"
    name=$(basename ${IID} _val_1.fq)
    echo ${name}
    cd ${VC_code}/2.Alignment_output/
    #Modify the job-name with the sample name
    R1_gz=${trimmed_dir}/${name}_val_1.fq
    R2_gz=${trimmed_dir}/${name}_val_2.fq
   sbatch --job-name="runAlign_${name}" --output="runAlign_${name}.out" --error="runAlign_${name}.err" ../scripts/alignment.V2.sh ${name} ${trimmed_dir} $R1_gz $R2_gz ${ref_genome} ${output_dir} ${my_softwares}
done

