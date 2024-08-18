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


out_dir=${ANGSD_results}/PCA
mkdir -p $out_dir


angsd_dir=${my_softwares}/angsd
angsd_out=${ANGSD_results}/ANGSD_out

pcangus_dir=/home/gabriel.zayas/.conda/envs/PCAngsd

mkdir -p ${ANGSD_code}/pca_out/
cd ${ANGSD_code}/pca_out/
for pop_info in $(ls ${angsd_out}/*.pop.info); do
    name=$(basename ${pop_info} .pop.info)
    beagle_file=${angsd_out}/${name}/${name}.beagle.gz
    N_ind=$(wc -l ${pop_info} | awk '{print $1}')
    echo "Processing $name with ${N_ind} in $beagle_file"
    echo ${name}
    sbatch --job-name="PCA_${name}" --output="PCA_${name}.out" --error="PCA_${name}.err" ../scripts/pca.sh ${ANGSD_code} ${pcangus_dir} ${name} ${pop_info} ${beagle_file} ${out_dir}/${name}
done


