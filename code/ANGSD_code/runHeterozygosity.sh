#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=Het_ANGSD_%j.out
#SBATCH --error=Het_ANGSD_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb


source ANGSD_project_env.sh
source ~/.bashrc
# Set the correct environment variables
export PATH=${my_softwares}/bin:$PATH
export LD_LIBRARY_PATH=${my_softwares}/lib:$LD_LIBRARY_PATH
export CFLAGS="-I${my_softwares}/include"
export LDFLAGS="-L${my_softwares}/lib"
# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $ANGSD_code"
echo "Data Directory: $ANGSD_data"
echo "Results Directory: $ANGSD_results"

# Create output directory
out_dir=${ANGSD_results}/Het/
mkdir -p $out_dir
rm $out_dir/*
# Define file paths
ref_genome=${ANGSD_data}/final_assembly/Hesperapis_oraria_2.curatedScaff.clean_sort_rename_nuc.fasta
angsd_dir=${my_softwares}/angsd
bam_dir=${proj_dir}/results/Variant_Calling/3.Cleaning/
bam_filelist=${out_dir}/bam.filelist
angsd_beagle_mafs=${ANGSD_results}/ANGSD_out/combined/combined.mafs.gz
mkdir ${ANGSD_code}/Heterozygosity_output/
# List BAM files
ls ${bam_dir}/*.marked_duplicates.bam > ${bam_filelist}
zcat ${angsd_beagle_mafs} | cut -f 1,2 | tail -n +2 > ${out_dir}/sites.txt
$angsd_dir/angsd sites index ${out_dir}/sites.txt
# Verify BAM file listing
if [ ! -s ${bam_filelist} ]; then
  echo "No BAM files found in ${bam_dir}"
  exit 1
fi

# Extract names of females
grep "female" ${meta_data} | awk -F "," '{print $1}' > ${out_dir}/Females.names

# Verify female names listing
if [ ! -s ${out_dir}/Females.names ]; then
  echo "No female names found in ${meta_data}"
  exit 1
fi

# Generate file list for females
grep -f ${out_dir}/Females.names ${bam_filelist} > ${out_dir}/Females.bam.filelist

# Verify female BAM file listing
if [ ! -s ${out_dir}/Females.bam.filelist ]; then
  echo "No BAM files found for females in ${bam_filelist}"
  exit 1
fi
cd ${ANGSD_code}/Heterozygosity_output/
# Estimate genome-wide heterozygosity
while read -r individual; do
  bam_file=$(grep "$individual" ${bam_filelist})
  
  if [ -z "$bam_file" ]; then
    echo "No BAM file found for $individual"
    continue
  fi

  echo "Processing $individual with BAM file $bam_file"
     sbatch --job-name="Hetero_$individual" --output="Hetero_$individual.out" --error="Hetero_$individual.err" ../scripts/heterozygosity.sh ${ANGSD_code} ${individual} ${bam_file} ${ref_genome} ${out_dir} ${out_dir}/sites.txt 
    
done < ${out_dir}/Females.names

echo "Genome-wide heterozygosity estimation completed."
