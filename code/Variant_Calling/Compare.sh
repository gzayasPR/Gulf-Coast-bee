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

deepvariant_vcf=${VC_data}/variants/DeepVariant/Hesperapis_oraria.vcf
gatk_vcf=${VC_data}/variants/GATK/Hesperapis_oraria.vcf
out_dir=${VC_results}/6.0.Compare_Calling/

mkdir -p $out_dir
ml bcftools
bcftools query -f '%ID \n' ${gatk_vcf} > ${out_dir}/gatk.SNPs.ID

bcftools query -f '%ID \n' ${deepvariant_vcf} > ${out_dir}/deep.SNPs.ID
comm -12 <(sort ${out_dir}/gatk.SNPs.ID) <(sort ${out_dir}/deep.SNPs.ID) > ${out_dir}/common.SNPs.ID
ml vcftools
vcftools --vcf ${gatk_vcf} --mac 1 --recode --out ${out_dir}/temp.GATK

bcftools query -f '%ID \n' ${out_dir}/temp.GATK.recode.vcf > ${out_dir}/gatk.temp.SNPs.ID

# Counting SNPs before filtering
common_SNP=$(wc -l ${out_dir}/common.SNPs.ID | awk '{print $1}')
gatk_SNP=$(wc -l ${out_dir}/gatk.temp.SNPs.ID | awk '{print $1}')
deep_SNP=$(wc -l ${out_dir}/deep.SNPs.ID | awk '{print $1}')

# Calculate the unique markers for each caller before filtering
gatk_only_SNP=$(($gatk_SNP - $common_SNP))
deep_only_SNP=$(($deep_SNP - $common_SNP))

# Calculate percentages before filtering
gatk_common_percent=$(echo "scale=4; ($common_SNP / $gatk_SNP) * 100" | bc)
gatk_only_percent=$(echo "scale=4; ($gatk_only_SNP / $gatk_SNP) * 100" | bc)
deep_common_percent=$(echo "scale=4; ($common_SNP / $deep_SNP) * 100" | bc)
deep_only_percent=$(echo "scale=4; ($deep_only_SNP / $deep_SNP) * 100" | bc)

# Output the results before filtering
echo "Markers in Common (before filtering): $common_SNP"
echo "Markers in GATK Only (before filtering): $gatk_only_SNP, which is $gatk_only_percent% of GATK markers"
echo "Markers in DeepVariant Only (before filtering): $deep_only_SNP, which is $deep_only_percent% of DeepVariant markers"
echo "Percentage of markers in common (before filtering):"
echo " - GATK: $gatk_common_percent%"
echo " - DeepVariant: $deep_common_percent%"
min_Q=20
min_meanDP=3
max_meanDP=55
max_locus_missing="0.75"
MAC=5
# Filter VCF file based on minimum depth
#vcftools --vcf ${out_dir}/temp.GATK.recode.vcf --max-missing ${max_locus_missing} --mac ${MAC} --min-meanDP ${min_meanDP} --minQ ${min_Q} --out ${out_dir}/fixed_missing.GATK --recode --recode-INFO-all

vcftools --vcf ${gatk_vcf} --max-missing ${max_locus_missing} --mac ${MAC} --min-meanDP ${min_meanDP} --max-meanDP ${max_meanDP} --minQ ${min_Q} --out ${out_dir}/fixed_missing.GATK --recode --recode-INFO-all
# Filter VCF file based on minimum depth
vcftools --vcf ${deepvariant_vcf} --max-missing ${max_locus_missing} --mac ${MAC} --min-meanDP ${min_meanDP} --max-meanDP ${max_meanDP} --minQ ${min_Q} --out ${out_dir}/fixed_missing.DeepVariant --recode --recode-INFO-all

bcftools query -f '%ID \n' ${out_dir}/fixed_missing.GATK.recode.vcf  > ${out_dir}/Filtered.gatk.temp.SNPs.ID
bcftools query -f '%ID \n' ${out_dir}/fixed_missing.DeepVariant.recode.vcf  > ${out_dir}/Filtered.deep.SNPs.ID
comm -12 <(sort ${out_dir}/Filtered.gatk.temp.SNPs.ID) <(sort ${out_dir}/Filtered.deep.SNPs.ID) > ${out_dir}/Filtered.common.SNPs.ID


Filtered_common_SNP=$(wc -l ${out_dir}/Filtered.common.SNPs.ID | awk '{print $1}')
Filtered_GATK_SNP=$(wc -l ${out_dir}/Filtered.gatk.temp.SNPs.ID | awk '{print $1}')
Filtered_DeepV_SNP=$( wc -l ${out_dir}/Filtered.deep.SNPs.ID| awk '{print $1}')

# Calculate the unique markers for each caller
GATK_only_SNP=$(($Filtered_GATK_SNP - $Filtered_common_SNP))
DeepV_only_SNP=$(($Filtered_DeepV_SNP - $Filtered_common_SNP))

# Calculate percentages
GATK_common_percent=$(echo "scale=4; ($Filtered_common_SNP / $Filtered_GATK_SNP) * 100" | bc)
GATK_only_percent=$(echo "scale=4; ($GATK_only_SNP / $Filtered_GATK_SNP) * 100" | bc)
DeepV_common_percent=$(echo "scale=4; ($Filtered_common_SNP / $Filtered_DeepV_SNP) * 100" | bc)
DeepV_only_percent=$(echo "scale=4; ($DeepV_only_SNP / $Filtered_DeepV_SNP) * 100" | bc)

# Output the results
echo "Markers in Common: $Filtered_common_SNP"
echo "Markers in GATK Only: $GATK_only_SNP, which is $GATK_only_percent% of GATK markers"
echo "Markers in DeepVariant Only: $DeepV_only_SNP, which is $DeepV_only_percent% of DeepVariant markers"
echo "Percentage of markers in common:"
echo " - GATK: $GATK_common_percent%"
echo " - DeepVariant: $DeepV_common_percent%"

VCF=${out_dir}/fixed_missing.GATK.recode.vcf 
cd $r_library
ml r/4.4.0
Rscript $proj_dir/code/Population_Genomics/scripts/PCA.V2.R ${out_dir} ${meta_data} $VCF 
VCF=${out_dir}/fixed_missing.DeepVariant.recode.vcf  
Rscript $proj_dir/code/Population_Genomics/scripts/PCA.V2.R ${out_dir} ${meta_data} $VCF 
