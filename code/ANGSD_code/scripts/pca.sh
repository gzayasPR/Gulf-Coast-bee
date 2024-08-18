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


ANGSD_code=$1
pcangus_dir=$2
name=$3
pop_info=$4
beagle_file=$5
out_dir=$6

cd ${ANGSD_code}

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

echo "Input 1 : Code Directory: $ANGSD_code"
echo "Input 2 : pca Code Directory: $pcangus_dir"
echo "Input 3 : Group Name: $name"
echo "Input 4 : Population Info: $pop_info"
echo "Input 5 : Beagle File: $beagle_file"
echo "Input 6 : Output_dir: $out_dir"


mkdir -p $out_dir

N_ind=$(wc -l ${pop_info} | awk '{print $1}')
af_den=$((${N_ind} * 2))
maf_num=2
maf=$(echo "scale=2; $maf_num / $af_den" | bc)
echo ${maf}
cd ${out_dir}
ml miniconda3
eval "$(conda shell.bash hook)"
conda activate ${pcangus_dir}
pcangsd -b ${beagle_file}  -t 8 --iter 1000 --maf ${maf} -o $out_dir/PCA --selection --sites_save

cd $r_library

ml r/4.4.0
Rscript $ANGSD_code/scripts/PCA.R ${out_dir} ${name} ${meta_data} ${pop_info}
Rscript $ANGSD_code/scripts/PCA.figure.R ${out_dir} ${name} ${meta_data} ${pop_info}
cd ${ANGSD_code}
