#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=Het_ANGSD_%j.out
#SBATCH --error=Het_ANGSD_%j.err
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

# Create output directory
out_dir=${ANGSD_results}/Het/

cd $r_library
ml r/4.4.0
Rscript $ANGSD_code/scripts/plot.Heterozygosity.R ${out_dir} ${meta_data}
cd ${ANGSD_code}
