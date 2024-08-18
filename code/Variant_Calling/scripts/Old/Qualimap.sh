#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --job-name=Quali_map
#SBATCH --output=Quali_map_%j.out
#SBATCH --error=Quali_map_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=120:00:00  # Adjusted walltime to 24 hours
#SBATCH --ntasks=48  # Using all available cores on a single node
#SBATCH --nodes=1  # Use a single node
#SBATCH --partition=atlas  # Selecting the atlas partition
#SBATCH --qos=normal  # Using the normal QoS
source VC_project_env.sh
# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $VC_code"
echo "Data Directory: $VC_data"
echo "Results Directory: $VC_results"
# Define the main output directory where all job-specific directories are located
output_dir=$VC_results/2.Alignment/
ml qualimap
cd ${output_dir}
echo "BLX2737	BLX2737.sorted_aligned_reads.bam	USA-AL_Baldwin_co	female
BLX2738	BLX2738.sorted_aligned_reads.bam	USA-AL_Baldwin_co	female
BLX2739	BLX2739.sorted_aligned_reads.bam	USA-AL_Baldwin_co	male
BLX2740	BLX2740.sorted_aligned_reads.bam	USA-AL_Baldwin_co	male
BLX2741	BLX2741.sorted_aligned_reads.bam	USA-AL_Baldwin_co	male
BLX2743	BLX2743.sorted_aligned_reads.bam	USA-FL_Escambia_co	female
BLX2744	BLX2744.sorted_aligned_reads.bam	USA-FL_Escambia_co	female
BLX2745	BLX2745.sorted_aligned_reads.bam	USA-FL_Escambia_co	male
BLX2747	BLX2747.sorted_aligned_reads.bam	USA-FL_Escambia_co	female
BLX2748	BLX2748.sorted_aligned_reads.bam	USA-FL_Escambia_co	female
BLX2751	BLX2751.sorted_aligned_reads.bam	USA-AL_Baldwin_co	male
BLX2752	BLX2752.sorted_aligned_reads.bam	USA-AL_Baldwin_co	male
BLX2753	BLX2753.sorted_aligned_reads.bam	USA-AL_Baldwin_co	female
BLX2754	BLX2754.sorted_aligned_reads.bam	USA-AL_Baldwin_co	female
BLX2755	BLX2755.sorted_aligned_reads.bam	USA-FL_Escambia_co	male
BLX2756	BLX2756.sorted_aligned_reads.bam	USA-FL_Escambia_co	male
BLX2757	BLX2757.sorted_aligned_reads.bam	USA-FL_Escambia_co	male
BLX2759	BLX2759.sorted_aligned_reads.bam	USA-FL_Escambia_co	female
" > all.info.txt
awk -F '\t' '{print $1,$2,$3}' all.info.txt > cleaned_site.location.txt
awk -F '\t' '{print $1,$2,$4}' all.info.txt > cleaned_sex.txt
awk -F '\t' '{print $1,$2}' all.info.txt > cleaned.txt
export JAVA_OPTS="-Xms1G -Xmx8G"
ml miniconda3 
eval "$(conda shell.bash hook)"

# Check if trim_galore is available
if ! command -v qualimap &> /dev/null; then
    echo "trim_galore could not be found in the conda environment ${VC_conda}"
    exit 1
fi
conda activate /project/beenome100/conda_envs/qualimap
qualimap multi-bamqc -r -d cleaned_site.location.txt -outfile site_multibamqc.pdf
qualimap multi-bamqc -r -d  cleaned_sex.txt -outfile sex_multibamqc.pdf
qualimap multi-bamqc -r -d  cleaned.txt -outfile multibamqc.pdf
conda deactivate

    