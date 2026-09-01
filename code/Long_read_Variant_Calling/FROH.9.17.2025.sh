#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --job-name=froh_analysis
#SBATCH --output=froh_analysis.out
#SBATCH --error=froh_analysis.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=8
#SBATCH --nodes=8
#SBATCH --mem-per-cpu=8gb

# Load environment and modules
source LR_VC_project_env.sh
ml vcftools
ml samtools
ml bcftools
ml miniconda3
plink_dir=${my_softwares}/plinkv1.9
eval "$(conda shell.bash hook)"
source ~/.bashrc

# Set directories and constants
RESULT_DIR_BASE="${LR_VC_results}/FROH"
MAX_ID_LENGTH=100 
MIN_Q=20
MIN_MQ=20
MIN_MEAN_DP=30
MAX_MEAN_DP=100
MAX_LOCUS_MISSING="0.95"
GQ_THRESHOLD=20
MIN_DP=30
MAX_DP=100
MIN_MQ=30  # Minimum Mapping Quality threshold
# Calculate genome length
GENOME_LENGTH=$(awk '{sum+=$2} END {print sum}' $ref_genome.fai)
echo "Total genome length calculated: $GENOME_LENGTH bp"

# Function to process a given VCF file and output directory
process_vcf_to_froh() {
    local output_dir=$1
    local vcf_file=$2
    mkdir -p ${output_dir}
    
    echo "Processing VCF file: ${vcf_file}"
    echo "Output directory: ${output_dir}"

    # Normalize and annotate VCF
    bcftools view $vcf_file | bcftools norm -m -any -f ${ref_genome} -Oz -o ${output_dir}/norm.vcf
    bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' ${output_dir}/norm.vcf | \
    awk -v max_len=$MAX_ID_LENGTH 'BEGIN{FS=OFS="\t"} 
    {
        if($3 != "." && length($3) > max_len) 
            $3 = substr($3, 1, max_len); 
        print
    }' >${output_dir}/norm_ann.vcf

    # Manually fix the header if needed (header_file.txt should be checked manually)
    bcftools view -h ${output_dir}/norm.vcf > ${output_dir}/header_file.txt
    bcftools reheader -h ${output_dir}/header_file.txt -o ${output_dir}/fixed_norm_ann.vcf ${output_dir}/norm_ann.vcf

    # Filter the VCF file by mapping quality and other criteria
    vcftools --vcf  ${output_dir}/norm_ann.vcf  --min-meanDP $MIN_MEAN_DP --max-meanDP $MAX_MEAN_DP \
             --minQ $MIN_Q --out ${output_dir}/norm_ann.filtered --recode --recode-INFO-all
    conda activate "${my_softwares}/plinkv1.9/env"
    # Convert VCF to PLINK format
    plink --vcf ${output_dir}/norm_ann.filtered.recode.vcf --make-bed --allow-extra-chr --out ${output_dir}/bee_hesperapis_plink

    # Run ROH analysis
    plink --bfile ${output_dir}/bee_hesperapis_plink --allow-extra-chr  --homozyg \
          --homozyg-snp 50 \
          --homozyg-kb 100 \
          --homozyg-density 50 \
          --homozyg-gap 500 \
          --homozyg-window-snp 50 \
          --homozyg-window-het 1 \
          --homozyg-window-missing 5 \
          --homozyg-window-threshold 0.05 \
          --allow-extra-chr \
          --out ${output_dir}/bee_hesperapis_roh

    # Calculate FROH for each individual
    awk -v genome_length=$GENOME_LENGTH '
    BEGIN{FS=OFS="\t"} 
    NR>1 {roh_length[$1]+=$8} 
    END {
        for (ind in roh_length) 
        print ind, roh_length[ind]/genome_length 
    }' ${output_dir}/bee_hesperapis_roh.hom > ${output_dir}/FROH_results.txt
    # After the --homozyg plink call, add:
    awk 'BEGIN{OFS="\t"} {print $1, $2}' ${ref_genome}.fai > ${output_dir}/chr_sizes.txt
    echo "FROH analysis completed. Results saved in ${output_dir}/FROH_results.txt"
}

# Process datasets for Female and Male samples
process_vcf_to_froh "${RESULT_DIR_BASE}/Female" "${LR_VC_results}/Variant.Calling/DeepVariant/Hesperapis_oraria.deepvariant.vcf.gz"
process_vcf_to_froh "${RESULT_DIR_BASE}/Male" "${LR_VC_results}/Variant.Calling/DeepVariant/Hesperapis_oraria_2.deepvariant.vcf.gz"