#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=Ne_ANGSD_%j.out
#SBATCH --error=Ne_ANGSD_%j.err
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


out_dir=${ANGSD_results}/Ne/
mkdir -p $out_dir
ref_genome=${ANGSD_data}/final_assembly/Hesperapis_oraria_2.curatedScaff.clean_sort_rename_nuc.fasta
angsd_dir=${my_softwares}/angsd
bam_dir=${proj_dir}/results/Variant_Calling/3.Cleaning/
bam_filelist=${out_dir}/bam.filelist

ls ${bam_dir}/*.marked_duplicates.bam  > ${bam_filelist}
cat ${bam_filelist}

cat ${bam_filelist} | xargs -n1 basename | awk -F "." '{print $1}' | sort | uniq > ${out_dir}/pop.info
pop_info=${out_dir}/pop.info
N_ind=$(wc -l ${pop_info} | awk '{print $1}')
# Define other variables
min_dp=$((3 * $N_ind))
max_dp=$((50 * $N_ind))
echo $N_ind
echo $min_dp
echo $max_dp
# Step 1: Generate the SAF file
$angsd_dir/angsd -bam ${bam_filelist} \
    -doSaf 1 \
    -gl 1 \
    -anc ${ref_genome} \
    -doMajorMinor 1 \
    -doMaf 1 \
    -doCounts 1 \
    -minQ 20 \
    -minMapQ 30 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -setMinDepth ${min_dp} \
    -setMaxDepth ${max_dp} \
    -minInd $N_ind \
    -P 8 \
    -out ${out_dir}/output

# Step 2: Estimate the Site Frequency Spectrum (SFS)
$angsd_dir/misc/realSFS ${out_dir}/output.saf.idx > ${out_dir}/output.sfs

# Step 3: Calculate Theta (θ) using the SFS
$angsd_dir/misc/realSFS saf2theta ${out_dir}/output.saf.idx -sfs ${out_dir}/output.sfs -outname ${out_dir}/output

$angsd_dir/misc/realSFS ${out_dir}/output.saf.idx -fold 1 -maxIter 100 -P 4 > ${out_dir}/smallFolded.sfs

cd $r_library
ml r/4.4.0
Rscript $ANGSD_code/scripts/plot.sfs.R ${out_dir}

# Step 4: Calculate effective population size (Ne)
mu=1.e-10  # Mutation rate per generation per site
generation_time=1  # Assume 5 years per generation, adjust as necessary
gunzip -c ${out_dir}/output.thetas.gz | awk '{sum+=$3} END {print sum}'

$angsd_dir/misc/ThetaStat print out.thetas.idx 2>/dev/null |head

# Extract the overall theta (sum of thetas for all bins) from the theta output
theta=$(gunzip -c ${out_dir}/output.thetas.gz | awk '{sum+=$3} END {print sum}')

# Calculate Ne using the formula: Ne = theta / (4 * mu)
Ne=$(echo "scale=5; ${theta} / (4 * ${mu})" | bc)

# Optionally, calculate Ne for a specific generation time
Ne_per_gen=$(echo "scale=5; ${Ne} / ${generation_time}" | bc)

# Print the results
echo "Effective Population Size (Ne): ${Ne}"
echo "Ne per Generation: ${Ne_per_gen}"