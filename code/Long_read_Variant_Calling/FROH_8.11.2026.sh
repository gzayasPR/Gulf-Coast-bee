#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --job-name=froh_plink
#SBATCH --output=froh_plink.%j.out
#SBATCH --error=froh_plink.%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gzayas97@ufl.edu


############################################
#            USER / ENV SETTINGS           #
############################################

# Project env (paths like $ref_genome, $LR_VC_results, $my_softwares, etc.)
source LR_VC_project_env.sh

# Tools / threads
THREADS=8
ml bcftools
ml vcftools
ml samtools
ml miniconda3
eval "$(conda shell.bash hook)"
source ~/.bashrc
plink_dir="${my_softwares}/plinkv1.9"   # path containing plink binary
conda activate "${plink_dir}/env" 2>/dev/null || true

# Inputs
VCF_IN="${LR_VC_results}/Variant.Calling/DeepVariant/Hesperapis_oraria.deepvariant.vcf.gz"
REF_FA="$ref_genome"                         # reference FASTA
FAI="${REF_FA}.fai"                          # FASTA index
BAM_IN="${LR_VC_results}/Alignment/Hesperapis_oraria.bam"                                    # optional: BAM for depth QC ("" to skip)

# Output base
OUT_BASE="${LR_VC_results}/FROH"
SAMPLE_SET="Female"                          # label for this run
OUT_DIR="${OUT_BASE}/${SAMPLE_SET}"
mkdir -p "${OUT_DIR}"

# Optional scaffolds control
# TODO: point this at a file listing the 17 chromosome-level scaffold names so the
# F_ROH denominator equals their summed length (report as assembled male-reference length).
KEEP_SCAFFOLDS_LIST=""                       # file with scaffolds to KEEP (one per line) or "" to keep all
EXCLUDE_SCAFFOLDS_LIST=""                    # file with scaffolds to EXCLUDE (one per line) or "" for none
MIN_SCAFF_LEN_BP=0                           # exclude scaffolds shorter than this (0 = keep all)

# Genotype / site filtering (tuned for high-coverage WGS)
MIN_GQ=30           # set genotypes with GQ < MIN_GQ to missing
MIN_DP=10           # set genotypes with DP < MIN_DP to missing
SITE_MIN_QUAL=30    # drop sites with QUAL < SITE_MIN_QUAL
MAX_MISSING=0.05    # drop sites with > MAX_MISSING missing (i.e., require callrate >= 0.95)
MIN_MQ=40           # require INFO/MQ >= MIN_MQ when present; ignored if missing

# PLINK ROH params (Profile A default; adjust after density check)
HOM_KB=200                 # --homozyg-kb
HOM_SNP=30                # --homozyg-snp; calibrate to observed NSNP in roh.hom
HOM_DENSITY=50            # --homozyg-density; default, appropriate for variant-only VCF
HOM_GAP=1000            # --homozyg-gap (kb break)
HOM_WIN_SNP=50            # --homozyg-window-snp
HOM_WIN_HET=1              # --homozyg-window-het (0 for high-cov DV is reasonable)
HOM_WIN_MISS=1             # --homozyg-window-missing
HOM_WIN_THRESH=0.05        # --homozyg-window-threshold
HOM_HET=3                  # --homozyg-het (max heterozygous calls per ROH)

# ROH length bins (bp) for classifying runs (keep consistent with your R script)
BIN_SHORT_MAX=500000       # <0.5 Mb
BIN_INTER_MAX=1000000      # 0.5–1 Mb
# >1 Mb is long by complement

############################################
#          LOG VERSIONS / SETTINGS         #
############################################
{
  echo "Date: $(date)"
  echo "Host: $(hostname)"
  echo "bcftools: $(bcftools --version 2>/dev/null | head -1)"
  echo "vcftools: $(vcftools --version 2>&1 | head -1)"
  echo "samtools: $(samtools --version 2>/dev/null | head -1)"
  echo "plink: $(${plink_dir}/plink --version 2>/dev/null | tr '\n' ' ')"
  echo "Threads: ${THREADS}"
  echo "VCF_IN: ${VCF_IN}"
  echo "REF_FA: ${REF_FA}"
  echo "KEEP_SCAFFOLDS_LIST: ${KEEP_SCAFFOLDS_LIST:-NONE}"
  echo "EXCLUDE_SCAFFOLDS_LIST: ${EXCLUDE_SCAFFOLDS_LIST:-NONE}"
  echo "MIN_SCAFF_LEN_BP: ${MIN_SCAFF_LEN_BP}"
  echo "Filters: MIN_GQ=${MIN_GQ} MIN_DP=${MIN_DP} SITE_MIN_QUAL=${SITE_MIN_QUAL} MAX_MISSING=${MAX_MISSING} MIN_MQ=${MIN_MQ}"
  echo "PLINK ROH: kb=${HOM_KB} snp=${HOM_SNP} density=${HOM_DENSITY} gap=${HOM_GAP} win_snp=${HOM_WIN_SNP} win_het=${HOM_WIN_HET} win_miss=${HOM_WIN_MISS} win_thresh=${HOM_WIN_THRESH} het=${HOM_HET}"
} | tee "${OUT_DIR}/run_params.log"

############################################
#            HELPER FUNCTIONS              #
############################################

calc_genome_len_bp() {
  # Optionally filter by KEEP/EXCLUDE and MIN_SCAFF_LEN_BP to define "autosomal genome length" baseline
  local fai="$1"
  local keep="${2:-}"
  local exclude="${3:-}"
  local minlen="${4:-0}"

  if [[ -n "$keep" && -s "$keep" ]]; then
    awk 'NR==FNR{keep[$1]; next} ($1 in keep){print $1"\t"$2}' "$keep" "$fai"
  else
    awk -v m="$minlen" '{if($2>=m) print $1"\t"$2}' "$fai"
  fi | awk -v ex="${exclude:-}" '
    BEGIN{
      if(ex!=""){
        while((getline < ex)>0){excl[$1]=1}
        close(ex)
      }
    }
    { if(!( $1 in excl )) {sum+=$2} }
    END{ print sum+0 }'
}

make_region_arg_from_list() {
  # Turn a newline list of scaffolds into a comma-separated string for bcftools -r
  awk '{printf "%s%s", (NR>1?",":""), $1}'
}

############################################
#          DEPTH OUTLIER SCANNING          #
############################################

HIGH_DEPTH_SCAFF="${OUT_DIR}/high_depth_scaffolds.txt"
if [[ -n "${BAM_IN}" && -s "${BAM_IN}" ]]; then
  echo "[INFO] Scanning BAM for high-depth scaffolds…" | tee -a "${OUT_DIR}/run_params.log"
  samtools coverage -o "${OUT_DIR}/coverage.tsv" "${BAM_IN}"
  awk 'NR>1{cov[NR]=$7; name[NR]=$1; sum+=$7; sum2+=$7*$7}
       END{
         n=NR-1; m=sum/n; v=(sum2/n - m*m); sd=(v>0?sqrt(v):0); thr=m+4*sd;
         for(i=2;i<=NR;i++) if(cov[i]>thr) print name[i]"\t"cov[i]"\t"m"\t"sd"\t"thr
       }' "${OUT_DIR}/coverage.tsv" > "${HIGH_DEPTH_SCAFF}" || true
  [[ -s "${HIGH_DEPTH_SCAFF}" ]] && echo "[WARN] High-depth scaffolds flagged: $(wc -l < ${HIGH_DEPTH_SCAFF})" | tee -a "${OUT_DIR}/run_params.log"
else
  : > "${HIGH_DEPTH_SCAFF}"  # create empty file
fi

############################################
#        DEFINE AUTOSOMAL GENOME LEN       #
############################################

# Build final exclude list = user exclude + high-depth
FINAL_EXCLUDE="${OUT_DIR}/exclude_scaffolds.txt"
if [[ -s "${EXCLUDE_SCAFFOLDS_LIST}" ]]; then
  cat "${EXCLUDE_SCAFFOLDS_LIST}" > "${FINAL_EXCLUDE}"
else
  : > "${FINAL_EXCLUDE}"
fi
# Append high-depth scaffolds if any
if [[ -s "${HIGH_DEPTH_SCAFF}" ]]; then
  cut -f1 "${HIGH_DEPTH_SCAFF}" >> "${FINAL_EXCLUDE}"
fi
# De-duplicate
if [[ -s "${FINAL_EXCLUDE}" ]]; then
  sort -u "${FINAL_EXCLUDE}" -o "${FINAL_EXCLUDE}"
fi

GENOME_LENGTH_BP=$(calc_genome_len_bp "${FAI}" "${KEEP_SCAFFOLDS_LIST:-}" "${FINAL_EXCLUDE:-}" "${MIN_SCAFF_LEN_BP}")
echo "AUTOSOMAL_GENOME_LENGTH_BP: ${GENOME_LENGTH_BP}" | tee -a "${OUT_DIR}/run_params.log"

############################################
#           VCF CLEANING STEPS             #
############################################

echo "[INFO] Normalizing and filtering VCF…" | tee -a "${OUT_DIR}/run_params.log"

# 1) Normalize & keep biallelic SNPs
bcftools norm -m -both -f "${REF_FA}" -Oz -o "${OUT_DIR}/norm.vcf.gz" "${VCF_IN}"
bcftools view -m2 -M2 -v snps -Oz -o "${OUT_DIR}/snps.vcf.gz" "${OUT_DIR}/norm.vcf.gz"
bcftools index -f "${OUT_DIR}/snps.vcf.gz"

# 2) Per-genotype cleaning: set low-quality genotypes to missing
bcftools +setGT "${OUT_DIR}/snps.vcf.gz" -Ou  -- \
  -t q -n . -i "FMT/GQ<${MIN_GQ} || FMT/DP<${MIN_DP} || FMT/DP=\".\"" \
| bcftools view -Oz -o "${OUT_DIR}/gtclean.vcf.gz" 
bcftools index -f "${OUT_DIR}/gtclean.vcf.gz"


# 3) Per-site filters: QUAL, callrate, MQ when present
# Detect if INFO/MQ is defined in the VCF header
if bcftools view -h "${OUT_DIR}/gtclean.vcf.gz" | grep -q 'ID=MQ,'; then
  # MQ present: apply QUAL, callrate, and MQ
  bcftools view  \
    -i "QUAL>=${SITE_MIN_QUAL} && F_MISSING<=${MAX_MISSING} && INFO/MQ>=${MIN_MQ}" \
    -Oz -o "${OUT_DIR}/sitefilt.vcf.gz" "${OUT_DIR}/gtclean.vcf.gz"
else
  # MQ absent: apply QUAL and callrate only
  bcftools view  \
    -i "QUAL>=${SITE_MIN_QUAL} && F_MISSING<=${MAX_MISSING}" \
    -Oz -o "${OUT_DIR}/sitefilt.vcf.gz" "${OUT_DIR}/gtclean.vcf.gz"
fi

bcftools index -f "${OUT_DIR}/sitefilt.vcf.gz"


# 4) Keep/exclude scaffolds and min scaffold length for the ROH call
VCF_FOR_ROH="${OUT_DIR}/for_roh.vcf.gz"
if [[ -n "${KEEP_SCAFFOLDS_LIST}" && -s "${KEEP_SCAFFOLDS_LIST}" ]]; then
  REGIONS=$(make_region_arg_from_list < "${KEEP_SCAFFOLDS_LIST}")
  bcftools view --threads "${THREADS}" -r "${REGIONS}" -Oz -o "${VCF_FOR_ROH}" "${OUT_DIR}/sitefilt.vcf.gz"
elif [[ -s "${FINAL_EXCLUDE}" || ${MIN_SCAFF_LEN_BP} -gt 0 ]]; then
  # Build a keep list from FAI minus excluded and below-length scaffolds
  awk -v m="${MIN_SCAFF_LEN_BP}" 'NR==FNR{ex[$1]=1; next} {if(!($1 in ex) && $2>=m) print $1}' \
      "${FINAL_EXCLUDE}" "${FAI}" > "${OUT_DIR}/keep_scaffolds.final.txt"
  REGIONS=$(make_region_arg_from_list < "${OUT_DIR}/keep_scaffolds.final.txt")
  bcftools view --threads "${THREADS}" -r "${REGIONS}" -Oz -o "${VCF_FOR_ROH}" "${OUT_DIR}/sitefilt.vcf.gz"
else
  ln -sf "sitefilt.vcf.gz" "${VCF_FOR_ROH}"
fi
bcftools index -f "${VCF_FOR_ROH}"

############################################
#        (Optional) SNP DENSITY CHECK      #
############################################

# Quick global SNP/Mb estimate for logging
TOTAL_SNPS=$(bcftools index -s "${VCF_FOR_ROH}" | awk 'NR>1{sum+=$3} END{print sum+0}')
if [[ -z "${TOTAL_SNPS}" ]]; then TOTAL_SNPS=0; fi
if [[ "${GENOME_LENGTH_BP}" -gt 0 ]]; then
  DENSITY_PER_MB=$(awk -v n="${TOTAL_SNPS}" -v g="${GENOME_LENGTH_BP}" 'BEGIN{if(g>0) printf "%.2f", (n/g)*1e6; else print 0}')
else
  DENSITY_PER_MB="NA"
fi
echo "Approx SNP density: ${DENSITY_PER_MB} SNPs/Mb over kept scaffolds" | tee -a "${OUT_DIR}/run_params.log"

############################################
#            PLINK CONVERSION              #
############################################

echo "[INFO] Converting to PLINK format…" | tee -a "${OUT_DIR}/run_params.log"
conda activate "${my_softwares}/plinkv1.9/env"
plink --vcf "${VCF_FOR_ROH}" --make-bed --allow-extra-chr --double-id \
  --threads "${THREADS}" --out "${OUT_DIR}/plink_data"

############################################
#                PLINK ROH                #
############################################

echo "[INFO] Running PLINK ROH…" | tee -a "${OUT_DIR}/run_params.log"
plink --bfile "${OUT_DIR}/plink_data" --allow-extra-chr --threads "${THREADS}" --homozyg \
  --homozyg-kb "${HOM_KB}" \
  --homozyg-snp "${HOM_SNP}" \
  --homozyg-density "${HOM_DENSITY}" \
  --homozyg-gap "${HOM_GAP}" \
  --homozyg-window-snp "${HOM_WIN_SNP}" \
  --homozyg-window-het "${HOM_WIN_HET}" \
  --homozyg-window-missing "${HOM_WIN_MISS}" \
  --homozyg-window-threshold "${HOM_WIN_THRESH}" \
  --homozyg-het "${HOM_HET}" \
  --out "${OUT_DIR}/roh"

# PLINK .hom columns: POS1 is field 7, POS2 is field 8 (bp); length_bp = POS2 - POS1 + 1
# Compute FROH (sum of ROH bp / autosomal genome length)
echo "[INFO] Computing FROH…" | tee -a "${OUT_DIR}/run_params.log"
awk -v gl="${GENOME_LENGTH_BP}" 'BEGIN{FS=OFS=" "}
     NR==1{next}
     {
       len_bp = ($8 - $7 + 1);   # POS1=field 7, POS2=field 8 in PLINK .hom
       s[$2]+=len_bp;   # $2 = IID (with --double-id, FID=IID by default)
     }
     END{
       if(gl<=0){gl=1}
       for (id in s) printf "%s\t%.8f\n", id, s[id]/gl
     }' "${OUT_DIR}/roh.hom" > "${OUT_DIR}/FROH_results.txt"
echo "FROH written: ${OUT_DIR}/FROH_results.txt" | tee -a "${OUT_DIR}/run_params.log"

############################################
#         ROH LENGTH CLASS SUMMARY         #
############################################

# Summarize counts/total length per class (short/intermediate/long) per individual
awk -v smax="${BIN_SHORT_MAX}" -v imax="${BIN_INTER_MAX}" 'BEGIN{FS=OFS=" "}
     NR==1{next}
     {
       len_bp = ($8 - $7 + 1);   # POS1=field 7, POS2=field 8 in PLINK .hom
       id=$2;
       if(len_bp < smax){ cls="short" }
       else if(len_bp <= imax){ cls="intermediate" }
       else { cls="long" }
       n[id,cls]++; L[id,cls]+=len_bp; total[id]+=len_bp;
     }
     END{
       printf "ID\tclass\tsegments\tlength_bp\n";
       for (k in n){
         split(k,a,SUBSEP); id=a[1]; cls=a[2];
         printf "%s\t%s\t%d\t%d\n", id, cls, n[k], L[id,cls];
       }
     }' "${OUT_DIR}/roh.hom" > "${OUT_DIR}/ROH_length_classes.txt"

echo "ROH class summary: ${OUT_DIR}/ROH_length_classes.txt" | tee -a "${OUT_DIR}/run_params.log"

echo "[DONE] $(date)"