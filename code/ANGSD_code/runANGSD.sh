#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=run_ANGSD_%j.out
#SBATCH --error=run_ANGSD_%j.err
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
echo "SoftwareDirectory: $my_softwares"
out_dir=${ANGSD_results}/ANGSD_out/
mkdir -p $out_dir

angsd_dir=${my_softwares}/angsd
bam_dir=${proj_dir}/results/Variant_Calling/3.Cleaning/
ref_genome=${ANGSD_data}/final_assembly/Hesperapis_oraria_2.curatedScaff.clean_sort_rename_nuc.fasta
combined_bam_filelist=${out_dir}/combined_bam.filelist

ls ${bam_dir}/*.marked_duplicates.bam  > ${combined_bam_filelist}
cat ${combined_bam_filelist}

# Separate BAM files based on sex
female_bam_filelist="${out_dir}/female_bam.filelist"
male_bam_filelist="${out_dir}/male_bam.filelist"

grep -w "female" ${meta_data} | awk -F "," '{print $1}' > ${out_dir}/Females.names
grep -w "male" ${meta_data} | awk -F "," '{print $1}' > ${out_dir}/Males.names

grep -f ${out_dir}/Females.names ${combined_bam_filelist} > ${female_bam_filelist}
grep -f ${out_dir}/Males.names ${combined_bam_filelist} > ${male_bam_filelist}

# Verify female BAM file listing
if [ ! -s ${female_bam_filelist} ]; then
  echo "No BAM files found for females"
  exit 1
fi

# Verify male BAM file listing
if [ ! -s ${male_bam_filelist} ]; then
  echo "No BAM files found for males"
  exit 1
fi

cd ${out_dir}
cat ${combined_bam_filelist} | xargs -n1 basename | awk -F "." '{print $1}' | sort | uniq > ${out_dir}/combined.pop.info
cat ${female_bam_filelist} | xargs -n1 basename | awk -F "." '{print $1}' | sort | uniq > ${out_dir}/female.pop.info
cat ${male_bam_filelist} | xargs -n1 basename | awk -F "." '{print $1}' | sort | uniq > ${out_dir}/male.pop.info

mkdir -p ${ANGSD_code}/ANGSD_output/

cd ${ANGSD_code}/ANGSD_output/

for bam_filelist in $(ls ${out_dir}/*_bam.filelist); do
    name=$(basename ${bam_filelist} _bam.filelist)
    N_ind=$(wc -l ${bam_filelist} | awk '{print $1}')
    echo "Processing $name with ${N_ind} in BAM file $bam_filelist"
    echo ${out_dir}
    sbatch --job-name="ANGSD_${name}" --output="ANGSD_${name}.out" --error="ANGSD_${name}.err" ../scripts/ANGSD.sh ${ANGSD_code} ${name} ${bam_filelist} ${N_ind} ${out_dir}/${name}
done
