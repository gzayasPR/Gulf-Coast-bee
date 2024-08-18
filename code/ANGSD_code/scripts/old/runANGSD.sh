#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=run_ANGSD_%j.out
#SBATCH --error=run_ANGSD_%j.err
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


out_dir=${ANGSD_results}/ANGSD_out
mkdir -p $out_dir

angsd_dir=${my_code}/angsd
bam_dir=${proj_dir}/results/Variant_Calling/3.Cleaning/
ref_genome=${ANGSD_data}/final_assembly/Hesperapis_oraria_2.curatedScaff.clean_sort_rename_nuc.fasta
bam_filelist=${out_dir}/bam.filelist


ls ${bam_dir}/*.marked_duplicates.bam  > ${bam_filelist}
cat ${bam_filelist}
cd ${out_dir}
cat ${bam_filelist} | xargs -n1 basename  | awk -F "." '{print $1}' | sort | uniq   > ${out_dir}/pop.info
#$angsd_dir/angsd -GL 1 -out genolike -nThreads 48 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam ${bam_filelist}
# Run ANGSD with filtering parameters
$angsd_dir/angsd  -bam ${bam_filelist} -GL 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -doMaf 1 -minQ 20 -minMapQ 30 -remove_bads 1 -only_proper_pairs 1 -setMinDepthInd 3 -setMaxDepthInd 50 -minMaf 0.01 -SNP_pval 1e-6 -P 8   -minInd 16  -out genolike