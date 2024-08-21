#!/bin/bash -l
#SBATCH --account=beenome100
#SBATCH --output=inbreeding_ANGSD_%j.out
#SBATCH --error=inbreeding_ANGSD_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=168:00:00
#SBATCH --partition=bigmem
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=48

source LR_VC_project_env.sh
source ~/.bashrc
# Set the correct environment variables
export PATH=${my_softwares}/bin:$PATH
export LD_LIBRARY_PATH=${my_softwares}/lib:$LD_LIBRARY_PATH
export CFLAGS="-I${my_softwares}/include"
export LDFLAGS="-L${my_softwares}/lib"
# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $LR_VC_code"
echo "Data Directory: $LR_VC_data"
echo "Results Directory: $LR_VC_results"


out_dir=${LR_VC_results}/inbreeding_ANGSD/
mkdir -p $out_dir
angsd_dir=${my_softwares}/angsd/env
ngsRelate_dir=${my_softwares}/ngsRelate
bam_dir=${LR_VC_results}/Alignment
bam_filelist=${out_dir}/bam.filelist
ls ${bam_dir}/*.bam  > ${bam_filelist}
cat ${bam_filelist}
N_ind=$(wc -l ${bam_filelist} | awk '{print $1}')

source ~/.bashrc
ml miniconda3
conda activate $angsd_dir
angsd -bam ${bam_filelist} \
    -dosnpstat 1 \
    -nind 2 \
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
    -setMinDepth 20 \
    -P 40  \
    -out ${out_dir}/combined
    
cat ${bam_filelist} | xargs -n1 basename  | awk -F "." '{print $1}'  > ${out_dir}/pop.info
### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat ${out_dir}/combined.mafs.gz | cut -f5 | sed 1d > ${out_dir}/freq

### run NgsRelate
$ngsRelate_dir/ngsRelate/ngsRelate -g ${out_dir}/combined.glf.gz -n 2 -f 1
cd ${LR_VC_code}