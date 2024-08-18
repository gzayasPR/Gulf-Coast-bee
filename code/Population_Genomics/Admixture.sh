#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=ADMIXTURE_%j.out
#SBATCH --error=ADMIXTURE_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb

source PG_project_env.sh
# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $PG_code"
echo "Data Directory: $PG_data"
echo "Results Directory: $PG_results"Í
out_dir=${PG_results}/ADMIXTURE/DeepVariant/
vcf_file=${PG_data}/variants/DeepVariant/Hesperapis_oraria_no.het.vcf.gz
meta_data=${PG_data}/variants/Samples.Metadata.csv
name=Hesperapis_oraria
min_dp="9"
min_locus_missing="0.15"
min_indiv_missing="0.15"
min_MAC="1"
prune="TRUE"
mkdir -p ${out_dir}
cd $r_library

ml r/4.4.0
Rscript $PG_code/scripts/Parameters.QC_LD.R ${out_dir} ${vcf_file} ${meta_data} ${name} ${min_dp} ${min_locus_missing} ${min_indiv_missing} ${min_MAC} ${prune}
cd ${out_dir}
ml plink2
plink2 --allow-extra-chr -vcf ${out_dir}/${name}_pruned.vcf.gz  --vcf-half-call h --make-bed --out ${out_dir}/admixture
cut -f 1 ${out_dir}/admixture.bim | sed -e 's/[^0-9]//g' > ${out_dir}/map2
cut -f 2- ${out_dir}/admixture.bim > ${out_dir}/map3
paste ${out_dir}/map2 ${out_dir}/map3 > ${out_dir}/admixture.bim

eval "$(conda shell.bash hook)"
conda activate /project/beenome100/conda_envs/admixture

for K in 1 2 3 4 5; \
do admixture --cv ${out_dir}/admixture.bed $K | tee log${K}.out; done

grep -h CV log*.out > ${out_dir}/CV.logs

cv_logs=${out_dir}/CV.logs
chosen_k=2
q_path=${out_dir}/admixture
ids=${out_dir}/admixture.fam

cd $r_library
Rscript $PG_code/scripts/plots.admixture.R ${out_dir} ${name} ${meta_data} ${cv_logs} ${chosen_k} ${q_path} ${ids}

################
out_dir=${PG_results}/ADMIXTURE/GATK/
vcf_file=${PG_data}/variants/GATK/Hesperapis_oraria.vcf
meta_data=${PG_data}/variants/Samples.Metadata.csv
name=Hesperapis_oraria
min_dp="9"
min_locus_missing="0.15"
min_indiv_missing="0.15"
min_MAC="1"
prune="TRUE"
mkdir -p ${out_dir}
cd $r_library

ml r/4.4.0
Rscript $PG_code/scripts/Parameters.QC_LD.R ${out_dir} ${vcf_file} ${meta_data} ${name} ${min_dp} ${min_locus_missing} ${min_indiv_missing} ${min_MAC} ${prune}
cd ${out_dir}
ml plink2
plink2 --allow-extra-chr -vcf ${out_dir}/${name}_pruned.vcf.gz  --vcf-half-call h --make-bed --out ${out_dir}/admixture
cut -f 1 ${out_dir}/admixture.bim | sed -e 's/[^0-9]//g' > ${out_dir}/map2
cut -f 2- ${out_dir}/admixture.bim > ${out_dir}/map3
paste ${out_dir}/map2 ${out_dir}/map3 > ${out_dir}/admixture.bim

eval "$(conda shell.bash hook)"
conda activate /project/beenome100/conda_envs/admixture

for K in 1 2 3 4 5; \
do admixture --cv ${out_dir}/admixture.bed $K | tee log${K}.out; done

grep -h CV log*.out > ${out_dir}/CV.logs

cv_logs=${out_dir}/CV.logs
chosen_k=2
q_path=${out_dir}/admixture
ids=${out_dir}/admixture.fam

cd $r_library
Rscript $PG_code/scripts/plots.admixture.R ${out_dir} ${name} ${meta_data} ${cv_logs} ${chosen_k} ${q_path} ${ids}

