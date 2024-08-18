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
name=$2
bam_filelist=$3
N_ind=$4
out_dir=$5

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
echo "Input 2 :  Name: ${name}"
echo "Input 3 : Bam file list: ${bam_filelist}"
echo "Input 4 : Number of Individuals: $N_ind"
echo "Input 5 : out_dir: ${out_dir}"


echo "Processing $name with ${N_ind} in BAM file $bam_filelist"
mkdir -p $out_dir
rm -f $out_dir/*
cp ${bam_filelist} $out_dir/bam.file
min_dp=$((3 * $N_ind))
max_dp=$((50 * $N_ind))

$angsd_dir/angsd  -bam $out_dir/bam.file \
    -dosnpstat 1 \
    -gl 1 \
    -doGlf 2 \
    -doHWE 1 \
    -domajorminor 1 \
    -snp_pval 1e-6 \
    -domaf 1 \
    -dogeno 3 \
    -dopost 2 \
    -hetbias_pval 0.1 \
    -minQ 20 \
    -minMapQ 30 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -doCounts 1 \
    -setMinDepth ${min_dp} \
    -setMaxDepth  ${max_dp} \
    -minMaf 0 \
    -minInd $N_ind \
    -P 8   \
    -out $out_dir/$name
