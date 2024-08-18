#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=inbreeding_%j.out
#SBATCH --error=inbreeding_%j.err
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


out_dir=${ANGSD_results}/inbreeding/
mkdir -p $out_dir
angsd_dir=${my_softwares}/angsd
ngsRelate_dir=${my_softwares}/ngsRelate
bam_dir=${proj_dir}/results/Variant_Calling/3.Cleaning/
cat ${bam_filelist}
bam_filelist="${ANGSD_results}/ANGSD_out/female/bam.file"
angsd_beagle_mafs=${ANGSD_results}/ANGSD_out/females/females.mafs.gz
N_ind=$(wc -l ${bam_filelist} | awk '{print $1}')
min_dp=$((3 * $N_ind))
max_dp=$((50 * $N_ind))
min_ind=$(($N_ind - 5))

$angsd_dir/angsd -bam ${bam_filelist} \
    -dosnpstat 1 \
    -doHWE 1 \
    -GL 2 \
    -doGlf 3 \
    -doCounts 1 \
    -doMajorMinor 1 \
    -doMaf 1 \
    -minQ 20 \
    -minMapQ 30 \
    -remove_bads 1 \
    -hetbias_pval 0.1 \
    -only_proper_pairs 1 \
    -setMinDepth ${min_dp} \
    -setMaxDepth ${max_dp} \
    -P 8   \
    -minInd ${min_ind} \
    -out ${out_dir}/combined
    
cat ${bam_filelist} | xargs -n1 basename  | awk -F "." '{print $1}'  > ${out_dir}/pop.info
### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat ${out_dir}/combined.mafs.gz | cut -f5 | sed 1d > ${out_dir}/freq

### run NgsRelate
$ngsRelate_dir/ngsRelate/ngsRelate -g ${out_dir}/combined.glf.gz -n 18 -f ${out_ds/plot.kinship.R ${out_dir} ${meta_data} ${ids}
cd ${ANGSD_code}