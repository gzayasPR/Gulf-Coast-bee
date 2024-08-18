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

source VC_project_env.sh

# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $VC_code"
echo "Data Directory: $VC_data"
echo "Results Directory: $VC_results"
echo "Conda Enviroment: $VC_conda"
short_reads_meta_data=${meta_data}
# Create a directory for Quality Control within results
mkdir -p ${VC_results}/1.Trim_QC

# Change to the Quality Control directory
cd ${VC_results}/1.Trim_QC
output_dir=${VC_results}/1.Trim_QC
mkdir -p "${output_dir}/pre_trim"
mkdir -p "${output_dir}/post_trim"
mkdir -p ${VC_code}/1.QC_output/

cd ${VC_code}/1.QC_output/
dos2unix $short_reads_meta_data

awk -F',' 'NR==1 {for (i=1; i<=NF; i++) {if ($i=="BioSample") b=i; if ($i=="fq1") f1=i; if ($i=="fq2") f2=i}} NR>1 {print $b","$f1","$f2}' "$short_reads_meta_data" > ${output_dir}/reads_columns.csv


# Read the CSV file line by line
while IFS=',' read -r BioSample fq1 fq2
do
    # Skip the header line
    if [ "$BioSample" != "BioSample" ]; then
        # Save the variables
        echo "Processing BioSample: $BioSample"
        echo "fq1: $fq1"
        echo "fq2: $fq2"
        R1=$fq1
        R2=$fq2
        cd ${VC_code}/1.QC_output/
        name=$BioSample
    # Modify the job-name with the sample name
    sbatch --job-name="runQC_${name}" --output="runQC_${name}.out" --error="runQC_${name}.err" ../scripts/trim_qualitycontrol.V2.sh ${name} $R1 $R2 ${output_dir} $my_softwares
    fi
done < "${output_dir}/reads_columns.csv"
