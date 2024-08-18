#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=Final_calling_%j.out
#SBATCH --error=Final_calling_%j.err
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

ref_genome=${VC_data}/final_assembly/Hesperapis_oraria_2.curatedScaff.clean_sort_rename_nuc.fasta

clean_dir=${VC_results}/4.Calling_V2/GATK/
output_dir=${VC_results}/4.5.Final_Calling_V2/
mkdir -p ${output_dir}
# Define the path to the gvcf_list.txt
gvcf_list=${output_dir}/gvcf_list.txt

# Remove existing gvcf_list.txt if it exists
rm -f $gvcf_list
ml bcftools
ml python
# Iterate over each GVCF file and add to the list
for IID in ${clean_dir}/*.g.vcf.gz; do
  name=$(basename ${IID} .g.vcf.gz)
   echo -e "${name}\t${IID}" >> $gvcf_list
   #bcftools index -f ${IID}
done
rm -f ${output_dir}/*.bed
# Generate the interval file (genome_intervals.bed or genome_intervals.list)
awk '{print $1"\t0\t"$2}' ${ref_genome}.fai > ${output_dir}/genome_intervals.bed

# Define the size threshold for each chunk (e.g., 55 MB)
size_threshold=55000000

# Function to split BED file into chunks based on size threshold
split_bed_by_size() {
    local input_bed=$1
    local output_prefix=$2
    local size_threshold=$3
    local current_size=0
    local chunk_index=0

    mkdir -p $(dirname $output_prefix)
    > ${output_prefix}_chunk_${chunk_index}.bed

    while read -r line; do
        chr=$(echo $line | awk '{print $1}')
        start=$(echo $line | awk '{print $2}')
        end=$(echo $line | awk '{print $3}')
        size=$((end - start))

        # Debugging: Print the values being read
#        printf "Processing: chr=%s, start=%d, end=%d, size=%d\n" "$chr" "$start" "$end" "$size"

        if (( current_size + size > size_threshold )); then
            chunk_index=$((chunk_index + 1))
            current_size=0
        fi

        echo -e "$line" >> ${output_prefix}_${chunk_index}.bed
        current_size=$((current_size + size))
    done < $input_bed
}


# Split the BED file based on cumulative size
split_bed_by_size ${output_dir}/genome_intervals.bed ${output_dir}/genome_intervals $size_threshold
rm ${output_dir}/genome_intervals_chunk_chunk_0.bed
rm ${output_dir}/genome_intervals_chunk_0.bed
mkdir -p ${VC_code}/genomicsdb_output

ml gatk

# Loop over each chunk and run GenomicsDBImport and GenotypeGVCFs in parallel
for chunk in ${output_dir}/genome_intervals_*.bed; do
    {
        name=$(basename $chunk .bed)
        echo $name
        chunk_output_dir=${output_dir}/genomicsdb_$(basename $chunk .bed)
        # Combine GVCF files using GenomicsDBImport
        cd ${VC_code}/genomicsdb_output 
        sbatch --job-name="chunk_${name}" --output="chunk_${name}.out" --error="chunk_${name}.err" ../scripts/joint_genotyping.sh  $name $chunk_output_dir $chunk $gvcf_list ${ref_genome} ${output_dir} ${VC_code}/VC_project_env.sh
    } 
done

