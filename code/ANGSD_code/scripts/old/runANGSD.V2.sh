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
out_dir=${ANGSD_results}/ANGSD_out
mkdir -p $out_dir

angsd_dir=${my_softwares}/angsd
bam_dir=${proj_dir}/results/Variant_Calling/3.Cleaning/
ref_genome=${ANGSD_data}/final_assembly/Hesperapis_oraria_2.curatedScaff.clean_sort_rename_nuc.fasta
bam_filelist=${out_dir}/bam.filelist

ls ${bam_dir}/*.marked_duplicates.bam  > ${bam_filelist}
cat ${bam_filelist}

# Separate BAM files based on sex
female_bam_filelist="${out_dir}/female_bam.filelist"
male_bam_filelist="${out_dir}/male_bam.filelist"

grep "female" ${meta_data} | awk -F "," '{print $1}' > ${out_dir}/Females.names
grep "male" ${meta_data} | awk -F "," '{print $1}' > ${out_dir}/Males.names

grep -f ${out_dir}/Females.names ${bam_filelist} > ${female_bam_filelist}
grep -f ${out_dir}/Males.names ${bam_filelist} > ${male_bam_filelist}

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
cat ${bam_filelist} | xargs -n1 basename  | awk -F "." '{print $1}' | sort | uniq > ${out_dir}/pop.info

# Run ANGSD for female diploids
echo "Processing females (diploids)"
$angsd_dir/angsd -bam "$female_bam_filelist" \
  -ref "$ref_genome" \
  -GL 1 \
  -doGlf 2 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -minQ 20 \
  -minMapQ 30 \
  -remove_bads 1 \
  -only_proper_pairs 1 \
  -setMinDepth 27 \
  -setMaxDepth 450 \
  -doCounts 1 \
  -minMaf 0 \
  -SNP_pval 1e-6 \
  -minInd 7 \
  -P 8 \
  -out ${out_dir}/females

if [ $? -ne 0 ]; then
  echo "ANGSD command failed for females"
  exit 1
fi

# Run ANGSD for male monoploids
echo "Processing males (monoploids)"
$angsd_dir/angsd -bam "$male_bam_filelist" \
  -GL 1 \
  -doGlf 2 \
  -doCounts 1 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -minQ 20 \
  -minMapQ 30 \
  -remove_bads 1 \
  -only_proper_pairs 1 \
  -setMinDepth 27 \
  -setMaxDepth 450 \
  -minMaf 0 \
  -SNP_pval 1e-6 \
  -minInd 7 \
  -P 8 \
  -out ${out_dir}/males

if [ $? -ne 0 ]; then
  echo "ANGSD command failed for males"
  exit 1
fi

$angsd_dir/angsd  -bam ${bam_filelist} -GL 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -doMaf 1 -minQ 20 -minMapQ 30 -remove_bads 1 -only_proper_pairs 1 -setMinDepth 54 -setMaxDepth 900  -minMaf 0.01 -SNP_pval 1e-6 -P 8   -minInd 16 -out ${out_dir}/combined