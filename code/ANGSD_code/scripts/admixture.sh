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


ANGSD_code=$1
name=$2
pop_info=$3
beagle_file=$4
N_ind=$5
reps=$6
ks=$7
chosen_k=$8
out_dir=$9

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
echo "Input 2 : Group Name: $name"
echo "Input 3 : Population Info: $pop_info"
echo "Input 4 : Beagle File: $beagle_file"
echo "Input 5 : Number of Samples: $N_ind"
echo "Input 6 : Number of reps: $reps"
echo "Input 7 : Number of ks: $ks"
echo "Input 8 : Chosen K: $chosen_k"
echo "Input 9: Output_dir: $out_dir"

ngs_dir=${my_softwares}/NGSadmix
mkdir -p $out_dir

source ~/.bashrc

af_den=$((${N_ind} * 2))
maf_num=2
maf=$(echo "scale=2; $maf_num / $af_den" | bc)
echo ${maf}
for K in $(seq 1 $ks); do
    for rep in $(seq 1 $reps); do
        $ngs_dir -likes ${beagle_file} -K $K -o ${out_dir}/myoutfiles.k$K.rep$rep -minMaf ${maf}
    done
done


cd $r_library

ml r/4.4.0
Rscript $ANGSD_code/scripts/plots.admixture.R ${out_dir} ${name} ${meta_data} ${chosen_k} ${pop_info} ${reps} ${ks}
Rscript $ANGSD_code/scripts/Figure.admixture.R ${out_dir} ${name} ${meta_data} ${chosen_k} ${pop_info} ${reps} ${ks}
cd ${ANGSD_code}

