#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=GD_%j.log
#SBATCH --error=GD_%j.log
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
out_dir=${PG_results}/Genetic_Diversit_plink4/${variant_caller}
vcf_file=${PG_data}/variants/${variant_caller}/Hesperapis_oraria.vcf
name=Hesperapis_oraria
king=${my_softwares}/king
ngsRelate_dir=${my_softwares}/ngsRelate

mkdir -p ${out_dir}
cd ${out_dir}
rm -f ${out_dir}/*

# NOTE: vcftools is intentionally NOT loaded. bcftools handles ploidy, sample
# subsetting, and site-level QUAL/depth filters. PLINK 1.9 handles missingness
# and the per-individual het/F calculation.
#
# KEY DIFFERENCE FROM THE PER-POPULATION VERSION:
#   All site-level filtering (biallelic, QUAL, depth, missingness) is now applied
#   ONCE to the joint 8-sample cohort, producing a SINGLE QC'd site set shared by
#   both populations. The FL/AL split happens only afterward, purely to compute
#   per-individual het/F. Consequences:
#     * ONE "markers remaining after QC" number (report this in Methods).
#     * Depth is averaged over all 8 samples, not per-population, so the low-
#       coverage female no longer pulls a small per-pop average across the window.
#     * N_SITES is identical for every female and N_MISS is 0 within each pop
#       (joint --geno drops any site missing in >=1 of 8), so denominators align.
#     * F remains per-population: --het runs on each --keep subset, so E(HOM)
#       allele frequencies are still estimated within-population (no Wahlund
#       inflation). Only the SITE SET is shared, not the frequencies.
ml bcftools
ml miniconda3
eval "$(conda shell.bash hook)"
# PLINK v1.9 lives in a conda env. The reshape of the .het output below depends
# on the v1.9 column layout (FID IID O(HOM) E(HOM) N(NM) F). plink2's --het
# columns differ (#FID IID O(HOM) E(HOM) OBS_CT F) and the awk would need editing.
plink_dir="${my_softwares}/plinkv1.9"
conda activate "${plink_dir}/env"

# ---- Filtering parameters ----
min_Q=20
min_meanDP=7
max_meanDP=10000
max_locus_missing=0.95     # keep sites genotyped in >= 95% of individuals
# PLINK --geno takes the max *missing* fraction, i.e. 1 - max_locus_missing:
geno_missing=$(awk -v m="${max_locus_missing}" 'BEGIN{printf "%g", 1 - m}')

# Monomorphic-site handling for the Ho denominator.
#   true  -> count sites monomorphic *within* a population in N_SITES, so
#            Ho = het / (all called sites). Matches vcftools --het (Ho LOWER).
#   false -> plink --het default: monomorphic sites dropped, so
#            Ho = het / (polymorphic sites only) (Ho HIGHER).
# F is mathematically identical either way; only Ho moves.
INCLUDE_MONOMORPHIC=true

# ---- Individuals to analyze ----
# Only these samples are extracted from the VCF; everything else is ignored.
KEEP_SAMPLES="BLX2738 BLX2747 BLX2753 BLX2754 BLX2743 BLX2744 BLX2748 BLX2759"

# ---- Force diploid genotypes (haploid males -> diploid) ----
bcftools +fixploidy ${vcf_file} -- -f 2 > ${out_dir}/diploid.vcf

# ---- Build keep list + per-population name lists (metadata INTERSECT keep-list) ----
printf "%s\n" ${KEEP_SAMPLES} | sort -u > ${out_dir}/keep.samples
grep "USA-FL_Escambia_co" ${meta_data} | awk -F "," '{print $1}' \
    | grep -Fxf ${out_dir}/keep.samples > ${out_dir}/FL.names
grep "USA-AL_Baldwin_co"  ${meta_data} | awk -F "," '{print $1}' \
    | grep -Fxf ${out_dir}/keep.samples > ${out_dir}/AL.names

# Sanity check: warn about any requested sample not matched to a population.
cat ${out_dir}/FL.names ${out_dir}/AL.names | sort -u > ${out_dir}/matched.samples
n_req=$(wc -l < ${out_dir}/keep.samples)
n_match=$(wc -l < ${out_dir}/matched.samples)
echo "[INFO] Requested ${n_req} samples; matched ${n_match} to FL/AL ($(wc -l < ${out_dir}/FL.names) FL, $(wc -l < ${out_dir}/AL.names) AL)."
if [[ "${n_req}" -ne "${n_match}" ]]; then
    echo "[WARN] These requested samples were NOT found in the metadata populations:"
    comm -23 ${out_dir}/keep.samples ${out_dir}/matched.samples | sed 's/^/  /'
fi

# =====================================================================
#  ALL SITE-LEVEL FILTERING APPLIED ONCE, ON THE JOINT COHORT
# =====================================================================

# ---- Subset to the analysis samples ONCE ----
bcftools view -S ${out_dir}/keep.samples \
    -Ov -o ${out_dir}/cohort.vcf ${out_dir}/diploid.vcf

# ---- biallelic + QUAL + depth, averaged over ALL kept samples ----
# AVG(FMT/DP) is now taken over the full cohort, not a per-population subset.
bcftools view -m2 -M2 \
    -e "QUAL<${min_Q} || AVG(FMT/DP)<${min_meanDP} || AVG(FMT/DP)>${max_meanDP}" \
    -Ov -o ${out_dir}/cohort_filtered.vcf ${out_dir}/cohort.vcf

# ---- Missingness filter ONCE -> THE definitive QC site set ----
# --geno over 8 samples drops any site missing in >=1 individual, so every
# retained site is genotyped in all 8 (=> zero within-pop missingness downstream).
plink --vcf ${out_dir}/cohort_filtered.vcf \
      --vcf-half-call m \
      --allow-extra-chr --double-id \
      --geno ${geno_missing} \
      --make-bed --out ${out_dir}/cohort_data

# ---- Per-stage attrition log (THE marker count to cite in Methods) ----
n_in=$(grep -vc '^#'   ${out_dir}/cohort.vcf)
n_bcf=$(grep -vc '^#'  ${out_dir}/cohort_filtered.vcf)
n_markers=$(wc -l < ${out_dir}/cohort_data.bim)
echo "[INFO] cohort input sites:                 ${n_in}"
echo "[INFO] post biallelic/QUAL/depth (bcftools): ${n_bcf}"
echo "[INFO] post-missingness (--geno, FINAL):     ${n_markers}"
echo "[INFO] ==> Markers remaining after QC (shared by FL & AL): ${n_markers}"

# =====================================================================
#  SPLIT PER POPULATION ONLY FOR PER-INDIVIDUAL het / F
# =====================================================================

# per-pop keep files (FID IID; --double-id made FID=IID)
awk '{print $1, $1}' ${out_dir}/FL.names > ${out_dir}/FL.keep
awk '{print $1, $1}' ${out_dir}/AL.names > ${out_dir}/AL.keep

for pop in FL AL; do
    # --het / --missing run on the SAME cohort_data variant set, restricted to
    # this population's individuals. --het re-estimates allele frequencies on the
    # kept subset, so E(HOM)/F stay within-population. Monomorphic-within-pop
    # sites are silently dropped by --het; the block below folds them back in.
    plink --bfile ${out_dir}/cohort_data --allow-extra-chr \
          --keep ${out_dir}/${pop}.keep \
          --het --out ${out_dir}/${pop}_site_heterozygosity
    plink --bfile ${out_dir}/cohort_data --allow-extra-chr \
          --keep ${out_dir}/${pop}.keep \
          --missing --out ${out_dir}/${pop}_site_missing

    if [[ "${INCLUDE_MONOMORPHIC}" == "true" ]]; then
        # Merge .het (poly-only: FID IID O(HOM) E(HOM) N(NM) F) with
        # .imiss (all sites: FID IID MISS_PHENO N_MISS N_GENO F_MISS) by IID.
        #   N_full = N_GENO - N_MISS      (called sites incl. monomorphic)
        #   het    = N(NM) - O(HOM)       (mono sites contribute 0 hets)
        #   O_full = N_full - het         (observed hom incl. monomorphic)
        #   E_full = E(HOM) + (N_full - N(NM))   (each mono site adds +1 exp hom)
        #   F      = plink F              (invariant to monomorphic sites)
        awk 'BEGIN{OFS="\t"}
             FNR==NR{ if(FNR==1) next;
                      Oh[$2]=$3; Eh[$2]=$4; Nnm[$2]=$5; Fv[$2]=$6; ord[++k]=$2; next }
             { if(FNR==1) next;
               nf=$5-$4; het=Nnm[$2]-Oh[$2]; of=nf-het; ef=Eh[$2]+(nf-Nnm[$2]);
               line[$2]=sprintf("%s\t%d\t%.5f\t%d\t%s",$2,of,ef,nf,Fv[$2]) }
             END{ print "INDV","O(HOM)","E(HOM)","N_SITES","F";
                  for(i=1;i<=k;i++) print line[ord[i]] }' \
            ${out_dir}/${pop}_site_heterozygosity.het \
            ${out_dir}/${pop}_site_missing.imiss \
            > ${out_dir}/${pop}_site_heterozygosity.vcftools.het
    else
        # plink --het default (monomorphic sites excluded from N_SITES).
        awk 'BEGIN{OFS="\t"}
             NR==1{print "INDV","O(HOM)","E(HOM)","N_SITES","F"; next}
             {print $2,$3,$4,$5,$6}' \
            ${out_dir}/${pop}_site_heterozygosity.het \
            > ${out_dir}/${pop}_site_heterozygosity.vcftools.het
    fi
done

# ---- Combine: AL (with header) then FL (header stripped), as in the original ----
(cat ${out_dir}/AL_site_heterozygosity.vcftools.het \
 && sed '1d' ${out_dir}/FL_site_heterozygosity.vcftools.het) \
 > ${out_dir}/${name}.het

# ---- Per-population sample lists ----
bcftools query -l ${out_dir}/cohort_data.fam >/dev/null 2>&1  # (fam is the record now)
cut -d' ' -f2 ${out_dir}/FL.keep > ${out_dir}/FL.pop_info
cut -d' ' -f2 ${out_dir}/AL.keep > ${out_dir}/AL.pop_info

# ---- Downstream R (unchanged; consumes the vcftools-format .het) ----
cd $r_library
ml r/4.4.3
Rscript $PG_code/scripts/Hetero.V3.R ${out_dir}/ ${out_dir}/${name}.het ${meta_data}