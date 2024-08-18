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

cd ${VC_results}/1.Trim_QC/
mkdir -p Check_Contamination
cd ${VC_results}/1.Trim_QC/Check_Contamination


ml python
ml conda
ml seqtk
ml blastn

conda create -n check_contamination python=3.9
conda activate check_contamination
conda install -c bioconda seqtk
conda install -c bioconda blast
conda install -c conda-forge biopython
conda install -c conda-forge matplotlib
# Create a directory for the BLAST database
mkdir -p ${VC_data}/blastdb

# Navigate to the directory
cd ${VC_data}/blastdb
# Download the 'nt' database
update_blastdb.pl --decompress nt

output_dir=${VC_results}/1.Trim_QC/Check_Contamination

for IID in ${trimmed_dir}*_R1.paired.fastq; do
    # Extract the base filename without the suffix "_R1.fastq.gz"
    name=$(basename ${IID} _R1.paired.fastq)
    echo ${name}
    cd ${VC_code}/2.Alignment_output/
    #Modify the job-name with the sampl
    python ${VC_code}/Filter_GC.py ${trimmed_dir}/${name}_R1.paired.fastq ${output_dir}/${name}_R1.paired_reads_gc_32.fastq ${output_dir}/${name}_R1.paired_reads_gc_41.fastq
        python ${VC_code}/Filter_GC.py ${trimmed_dir}/${name}_R2.paired.fastq ${output_dir}/${name}_R2.paired_reads_gc_32.fastq ${output_dir}/${name}_R2.paired_reads_gc_41.fastq
        seqtk seq -a reads_gc_32.fastq > reads_gc_32.fasta
seqtk seq -a reads_gc_41.fastq > reads_gc_41.fasta
blastn -query reads_gc_32.fasta -db nt -out results_gc_32.txt -outfmt 6 -max_target_seqs 10
blastn -query reads_gc_41.fasta -db nt -out results_gc_41.txt -outfmt 6 -max_target_seqs 10

done


        seqtk seq -a reads_gc_32.BLX2759_R2.paired.fastq > reads_gc_32.BLX2759_R2.paired.fasta
        seqtk seq -a reads_gc_41.BLX2759_R2.paired.fastq > reads_gc_41.BLX2759_R2.paired.fasta
blastn -query reads_gc_32.BLX2759_R2.paired.fasta -db nt -out results_gc_32.txt -outfmt 6 -max_target_seqs 10
blastn -query reads_gc_41.BLX2759_R2.paired.fasta -db nt -out results_gc_41.txt -outfmt 6 -max_target_seqs 10
