#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --job-name="$1"
#SBATCH --output=$1.out
#SBATCH --error=$1.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=72:00:00  # Adjusted walltime to 24 hours
#SBATCH --ntasks=48  # Using all available cores on a single node
#SBATCH --nodes=1  # Use a single node
#SBATCH --partition=atlas  # Selecting the atlas partition
#SBATCH --qos=normal  # Using the normal QoS

ANGSD_code=$1
individual=$2
bam_file=$3
ref_genome=$4
out_dir=$5
sites=$6
cd ${ANGSD_code}

source ANGSD_project_env.sh
source ~/.bashrc
# Set the correct environment variables
export PATH=${my_softwares}/bin:$PATH
export LD_LIBRARY_PATH=${my_softwares}/lib:$LD_LIBRARY_PATH
export CFLAGS="-I${my_softwares}/include"
export LDFLAGS="-L${my_softwares}/lib"
angsd_dir=${my_softwares}/angsd

echo "Input 1 : Code Directory: $ANGSD_code"
echo "Input 2 : Individual Name: $individual"
echo "Input 3 : Bam file: $bam_file"
echo "Input 4 : Ref genome: $ref_genome"
echo "Input 5 : Output_dir: $out_dir"

if [ -z "$bam_file" ]; then
   echo "No BAM file found for $individual"
   continue
fi
echo "Processing $individual with BAM file $bam_file"
rm "${out_dir}/$individual.saf.idx"
rm "${out_dir}/$individual.saf.gz"
rm "${out_dir}/$individual.saf.pos.gz"
rm "${out_dir}/$individual.mafs.gz"
# Run ANGSD to calculate SAF with quality control and call non-polymorphic sites
  $angsd_dir/angsd -i "$bam_file" \
    -anc "$ref_genome" \
    -ref "$ref_genome" \
    -dosaf 1 \
    -GL 1 \
    -P 8 \
    -C 50 \
    -minQ 20 \
    -minMapQ 30 \
    -sites $sites \
   -setMinDepth 7 \
   -setMaxDepth 50 \
   -doCounts 1 \
   -only_proper_pairs 1 \
   -remove_bads 1 \
    -out "${out_dir}/$individual"

if [ $? -ne 0 ]; then
   echo "ANGSD command failed for $individual"
    continue
  fi
rm "${out_dir}/$individual.est.ml"
# Estimate heterozygosity
$angsd_dir/misc/realSFS -fold 1 "${out_dir}/$individual.saf.idx"  > "${out_dir}/$individual.est.ml"

  if [ $? -ne 0 ]; then
    echo "realSFS command failed for $individual"
  fi
