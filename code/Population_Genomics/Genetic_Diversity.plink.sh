#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=GD_%j.out
#SBATCH --error=GD_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb

source PG_project_env.sh
source ~/.bashrc

# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $PG_code"
echo "Data Directory: $PG_data"
echo "Results Directory: $PG_results"

# Set the correct environment variables
export PATH=${my_softwares}/bin:$PATH
export LD_LIBRARY_PATH=${my_softwares}/lib:$LD_LIBRARY_PATH
export CFLAGS="-I${my_softwares}/include"
export LDFLAGS="-L${my_softwares}/lib"

variant_caller=DeepVariant
out_dir=${PG_results}/Genetic_Diversity_plink/${variant_caller}
vcf_file=${PG_data}/variants/${variant_caller}/Females_Hesperapis_oraria.vcf
# meta_data=${PG_data}/variants/Samples.Metadata.csv   # defined in PG_project_env.sh
name=Hesperapis_oraria_Hetero
king=${my_softwares}/king
ngsRelate_dir=${my_softwares}/ngsRelate

mkdir -p ${out_dir}
cd ${out_dir}
rm -f ${out_dir}/*

# NOTE: vcftools is intentionally NOT loaded anywhere in this pipeline.
# bcftools handles ploidy, sample subsetting, and site-level QUAL/depth filters.
# PLINK 1.9 handles the missingness filter and the per-individual het/F calculation.
ml bcftools
ml miniconda3
eval "$(conda shell.bash hook)"
# PLINK v1.9 lives in a conda env (same convention as FROH.sh). The reshape of
# the .het output below depends on the v1.9 column layout
# (FID IID O(HOM) E(HOM) N(NM) F). If you switch to plink2 (ml plink2), its
# --het columns differ (#FID IID O(HOM) E(HOM) OBS_CT F) and the awk must change.
plink_dir="${my_softwares}/plinkv1.9"   # path containing the plink v1.9 conda env
conda activate "${plink_dir}/env"

# ---- Filtering parameters (unchanged from the vcftools version) ----
min_Q=20
min_meanDP=3
max_meanDP=55
max_locus_missing=0.95      # keep sites genotyped in >= 95% of individuals
# PLINK --geno takes the max *missing* fraction, i.e. 1 - max_locus_missing:
geno_missing=$(awk -v m="${max_locus_missing}" 'BEGIN{printf "%g", 1 - m}')

# The following were declared in the original script but never applied to any
# filtering command, so they are left unused here to reproduce identical results.
# Uncomment / wire them in only if you actually want per-genotype DP/GQ filtering
# (note: minDP=10 would drop many calls on the ~15.8x-coverage individual).
# GQ_threshold=20
# minDP=10
# maxDP=100

# ---- Force diploid genotypes (haploid males -> diploid) ----
bcftools +fixploidy ${vcf_file} -- -f 2 > ${out_dir}/diploid.vcf

# ---- Build per-site sample lists from metadata (unchanged) ----
grep "USA-FL_Escambia_co" ${meta_data} | awk -F "," '{print $1}' > ${out_dir}/FL.names
grep "USA-AL_Baldwin_co"  ${meta_data} | awk -F "," '{print $1}' > ${out_dir}/AL.names

# ---- Split into site-specific VCFs (bcftools view -S replaces vcftools --keep) ----
bcftools view -S ${out_dir}/FL.names --force-samples -Ov -o ${out_dir}/FL_site.vcf ${out_dir}/diploid.vcf
bcftools view -S ${out_dir}/AL.names --force-samples -Ov -o ${out_dir}/AL_site.vcf ${out_dir}/diploid.vcf

# ---- Per-population site-level filtering (bcftools replaces vcftools) ----
# biallelic (-m2 -M2) + QUAL >= min_Q + mean depth within [min_meanDP, max_meanDP].
# AVG(FMT/DP) is bcftools' equivalent of vcftools --min-meanDP / --max-meanDP.
# Add -v snps to the view call if you want to restrict to SNPs (vcftools --het did not).
for pop in FL AL; do
    bcftools view -m2 -M2 \
        -e "QUAL<${min_Q} || AVG(FMT/DP)<${min_meanDP} || AVG(FMT/DP)>${max_meanDP}" \
        -Ov -o ${out_dir}/${pop}_site_filtered.vcf \
        ${out_dir}/${pop}_site.vcf
done

# ---- Per-individual heterozygosity / F (PLINK --het replaces vcftools --het) ----
# --geno applies the missingness threshold (replaces vcftools --max-missing).
# --allow-extra-chr keeps the non-integer bee contig/scaffold names (verified:
#   these ARE included in --het). --double-id keeps full sample IDs as FID=IID.
# Allele frequencies for E(HOM) are estimated within each population's file,
# matching the per-population behavior of the vcftools version.
for pop in FL AL; do
    plink --vcf ${out_dir}/${pop}_site_filtered.vcf \
            -vcf-half-call m \
          --allow-extra-chr \
          --double-id \
          --geno ${geno_missing} \
          --het \
          --out ${out_dir}/${pop}_site_heterozygosity

    # Reshape PLINK's .het (FID IID O(HOM) E(HOM) N(NM) F) into the vcftools
    # .het layout (INDV O(HOM) E(HOM) N_SITES F) that Hetero.V3.R expects.
    awk 'BEGIN{OFS="\t"}
         NR==1{print "INDV","O(HOM)","E(HOM)","N_SITES","F"; next}
         {print $2,$3,$4,$5,$6}' \
        ${out_dir}/${pop}_site_heterozygosity.het \
        > ${out_dir}/${pop}_site_heterozygosity.vcftools.het
done

# ---- Combine: AL (with header) then FL (header stripped), as in the original ----
(cat ${out_dir}/AL_site_heterozygosity.vcftools.het \
 && sed '1d' ${out_dir}/FL_site_heterozygosity.vcftools.het) \
 > ${out_dir}/${name}.het

# ---- Per-population sample lists (bcftools query -l, unchanged) ----
bcftools query -l ${out_dir}/FL_site_filtered.vcf > ${out_dir}/FL.pop_info
bcftools query -l ${out_dir}/AL_site_filtered.vcf > ${out_dir}/AL.pop_info

# ---- Downstream R (unchanged; consumes the vcftools-format .het) ----
cd $r_library
ml r/4.4.3
Rscript $PG_code/scripts/Hetero.V3.R ${out_dir}/ ${out_dir}/${name}.het ${meta_data}