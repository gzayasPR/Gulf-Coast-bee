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


out_dir=${ANGSD_results}/ADMIXTURE/Females/
mkdir -p $out_dir

angsd_dir=${my_softwares}/angsd
angsd_out=${ANGSD_results}/ANGSD_out/females.beagle.gz
pop_info=${ANGSD_results}/ANGSD_out/pop.info
bam_filelist=${ANGSD_results}/ANGSD_out/female_bam.filelist
cat ${bam_filelist} | xargs -n1 basename  | awk -F "." '{print $1}' | sort | uniq > ${out_dir}/pop.info

cd ${out_dir}

for K in 1 2 3 4 5; \
    do
    for rep in 1 2 3 4 5 6 7 8 9 10 ; \
        do $angsd_dir/misc/NGSadmix -likes ${angsd_out} -K $K -o ${out_dir}/myoutfiles.k$K.rep$rep -minMaf 0.05 -minInd 8
        done; done

name=Females_Hesperapis_oraria
chosen_k=2
reps=10
ks=5
ids=${out_dir}/pop.info
cd $r_library
ml r/4.4.0
Rscript $ANGSD_code/scripts/plots.admixture.R ${out_dir} ${name} ${meta_data} ${chosen_k} ${ids} ${reps} ${ks}


Rscript $ANGSD_code/scripts/Figure.admixture.R ${out_dir} ${name} ${meta_data} ${chosen_k} ${ids} ${reps} ${ks}
cd ${ANGSD_code}

