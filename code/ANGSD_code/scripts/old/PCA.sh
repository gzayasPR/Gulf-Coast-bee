#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=PCA_ANGSD_%j.out
#SBATCH --error=PCA_ANGSD_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb

source ANGSD_project_env.sh
# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $ANGSD_code"
echo "Data Directory: $ANGSD_data"
echo "Results Directory: $ANGSD_results"


out_dir=${ANGSD_results}/PCA
mkdir -p $out_dir

angsd_dir=${my_softwares}/angsd
pcangus_dir=/home/gabriel.zayas/.conda/envs/PCAngsd
angsd_beagle=${ANGSD_results}/ANGSD_out_V5/combined.beagle.gz
pop_info=${ANGSD_results}/ANGSD_out/pop.info

cd ${out_dir}
ml miniconda3
eval "$(conda shell.bash hook)"
conda activate ${pcangus_dir}
pcangsd -b ${angsd_beagle} -t 8 --iter 1000 --maf 0.125 -o $out_dir/PCA --selection --sites_save

cd $r_library


name=Hesperapis_oraria
ids=${pop_info}
cd $r_library
ml r/4.4.0
Rscript $ANGSD_code/scripts/PCA.R ${out_dir} ${name} ${meta_data} ${ids}
Rscript $ANGSD_code/scripts/PCA.figure.R ${out_dir} ${name} ${meta_data} ${ids}
cd ${ANGSD_code}
