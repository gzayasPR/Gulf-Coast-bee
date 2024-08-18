#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=3.5.Intial_calling_%j.out
#SBATCH --error=3.5.Intial_calling_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb

source VC_project_env.sh
# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $VC_code"
echo "Data Directory: $VC_data"
echo "Results Directory: $VC_results"

ref_genome=${VC_data}/final_assembly/Hesperapis_oraria_2.curatedScaff.clean_sort_rename_nuc.fasta

clean_dir=${VC_results}/3.Cleaning/
output_dir=${VC_results}/3.5.Intial_calling/
mkdir -p ${output_dir}
metrics=${output_dir}/metrics.txt
rm -f $metrics
ml gatk
# Iterate over each GVCF file and add to the list
for IID in ${clean_dir}/*.raw_variants.vcf; do
    name=$(basename ${IID} .g.vcf.gz)
gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=48" VariantFiltration \
    -R ${ref_genome} \
    -V ${IID} \
    -O ${output_dir}/${name}.initial.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "FAIL"
    vcf_total=$(wc -l ${output_dir}/${name}.initial.vcf)
    vcf_passes=$(grep "PASS" ${output_dir}/${name}.initial.vcf | wc-l)
    echo -e "${name}\t${vcf_total}\t${vcf_passes}" >> $metrics
done
