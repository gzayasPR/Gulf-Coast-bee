#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=ADMIXTURE_ANGSD_%j.out
#SBATCH --error=ADMIXTURE_ANGSD_%j.err
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


out_dir=${ANGSD_results}/ADMIXTURE
mkdir -p $out_dir


angsd_dir=${my_softwares}/angsd
angsd_out=${ANGSD_results}/ANGSD_out

chosen_k=2
reps=10
ks=5

mkdir -p ${ANGSD_code}/admixture_out/
cd ${ANGSD_code}/admixture_out/

for pop_info in $(ls ${angsd_out}/*pop.info); do
    name=$(basename ${pop_info} .pop.info)
    beagle_file=${angsd_out}/${name}/${name}.beagle.gz
    N_ind=$(wc -l ${pop_info} | awk '{print $1}')
    echo "Processing $name with ${N_ind} in $beagle_file"
    echo ${out_dir}
    sbatch --job-name="ADMIXTURE_${name}" --output="ADMIXTURE_${name}.out" --error="ADMIXTURE_${name}.err" ../scripts/admixture.sh ${ANGSD_code} ${name} ${pop_info} ${beagle_file} ${N_ind} ${reps} ${ks} ${chosen_k} ${out_dir}/${name}
done


