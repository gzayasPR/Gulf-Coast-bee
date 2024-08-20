#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=MergingVCF_%j.out
#SBATCH --error=MergingVCF_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb

source VC_project_env.sh
# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $VC_code"
echo "Data Directory: $VC_data"
echo "Results Directory: $VC_results"

cd $VC_code
GLNexus_dir=${my_softwares}/GLNexus/
vcf_dir=${VC_results}/4.Calling/DeepVariant/
output_dir=${VC_results}/5.Merge_Calling/DeepVariant
mkdir -p ${output_dir}
# Define the path to the gvcf_list.txt
g_vcf_list=${output_dir}/g.vcf_list.txt
 
rm -f $g_vcf_list
# Iterate over each GVCF file and add to the list

for vcf in $(ls ${vcf_dir}/*.deepvariant.g.vcf.gz); do
    echo -e "${vcf}" >> $g_vcf_list
done
ml apptainer
TMPDIR=$output_dir/temp
mkdir -p $TMPDIR
export APPTAINER_CACHEDIR=$TMPDIR 
export APPTAINER_TMPDIR=$TMPDIR
VERSION="1.6.1"
cd ${output_dir}
rm -rf ${output_dir}/GLnexus.DB

apptainer run \
  ${GLNexus_dir}/glnexus_v1.2.7.sif \
  /usr/local/bin/glnexus_cli \
  --config DeepVariantWGS \
  ${vcf_dir}/*.deepvariant.g.vcf.gz > ${output_dir}/deepvariant.cohort.bcf


ml bcftools  
bcftools convert -O v -o ${output_dir}/deepvariant.cohort.vcf ${output_dir}/deepvariant.cohort.bcf
bcftools view ${output_dir}/deepvariant.cohort.vcf| bcftools norm -m -any -f ${ref_genome} -Oz -o ${output_dir}/norm.deepvariant.cohort.vcf
mkdir -p ${VC_data}/variants/DeepVariant/
cp ${output_dir}/norm.deepvariant.cohort.vcf ${VC_data}/variants/DeepVariant/temp.vcf

bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' ${VC_data}/variants/DeepVariant/temp.vcf > ${VC_data}/variants/DeepVariant/Hesperapis_oraria.vcf

grep "female" $meta_data | awk -F "," '{print $1}' > ${VC_data}/variants/DeepVariant/Females.ID

bcftools view -S ${VC_data}/variants/DeepVariant/Females.ID ${VC_data}/variants/DeepVariant/Hesperapis_oraria.vcf > ${VC_data}/variants/DeepVariant/Females_Hesperapis_oraria.vcf
