#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=KIN_ANGSD_%j.out
#SBATCH --error=KIN_ANGSD_%j.err
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


out_dir=${ANGSD_results}/kinship/
mkdir -p $out_dir
angsd_dir=${my_softwares}/angsd
ngsRelate_dir=${my_softwares}/ngsRelate
bam_dir=${proj_dir}/results/Variant_Calling/3.Cleaning/
bam_filelist=${out_dir}/bam.filelist

ls ${bam_dir}/*.marked_duplicates.bam  > ${bam_filelist}
cat ${bam_filelist}
angsd_beagle_mafs=${ANGSD_results}/ANGSD_out/combined/combined.mafs.gz
zcat ${angsd_beagle_mafs} | cut -f 1,2 | tail -n +2 > ${out_dir}/sites.txt
$angsd_dir/angsd sites index ${out_dir}/sites.txt

$angsd_dir/angsd -bam ${bam_filelist} \
    -GL 2 \
    -doGlf 3 \
    -doCounts 1 \
    -doMajorMinor 1 \
    -doMaf 1 \
    -sites ${out_dir}/sites.txt \
    -minQ 20 \
    -minMapQ 30 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -setMinDepth 54 \
    -setMaxDepth 900  \
    -minMaf 0.125 \
    -SNP_pval 1e-6 \
    -P 8   \
    -minInd 16 \
    -out ${out_dir}/combined
    
cat ${bam_filelist} | xargs -n1 basename  | awk -F "." '{print $1}'  > ${out_dir}/pop.info
### Then we extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat ${out_dir}/combined.mafs.gz | cut -f5 |sed 1d > ${out_dir}/freq

### run NgsRelate
$ngsRelate_dir/ngsRelate/ngsRelate -g ${out_dir}/combined.glf.gz -n 18 -f ${out_dir}/freq  -O ${out_dir}/newres
$ngsRelate_dir/ngsRelate/ngsRelate -g ${out_dir}/combined.glf.gz -n 18 -f ${out_dir}/freq -F 1 -O ${out_dir}/inbreeding

ids=${out_dir}/pop.info
cd $r_library
ml r/4.4.0
Rscript $ANGSD_code/scripts/plot.kinship.R ${out_dir} ${meta_data} ${ids}
cd ${ANGSD_code}
