#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --job-name=psmc_hesp
#SBATCH --output=psmc_hesp_%j.log
#SBATCH --error=psmc_hesp_%j.log
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zayas@wisc.edu

# =============================================================================
# PSMC Historical Ne Pipeline  (single diploid female)  -- LOOPING VERSION
#
# Reconstructs historical effective population size from one high-coverage
# diploid genome using the Pairwise Sequentially Markovian Coalescent.
#   Li & Durbin 2011, Nature 475:493-496 ; GitHub: https://github.com/lh3/psmc
#
# WHAT CHANGED VS THE SINGLE-RUN VERSION:
#   MU and PSMC_PATTERN are now ARRAYS (MU_LIST, PSMC_PATTERN_LIST). The script
#   sweeps every (pattern, mu) combination in one submission, reusing work at
#   exactly the level each parameter touches:
#
#     Steps 1-4  (shared)      run ONCE, before any loop. Independent of
#                              MU / GEN / pattern.
#     Steps 5-6  (per pattern) OUTER loop over PSMC_PATTERN_LIST. The psmc
#                              binary + bootstrap define the MODEL, so each
#                              pattern gets its own models/<tag>/ built once.
#     Step 7     (per mu)      INNER loop over MU_LIST. MU is PLOT-ONLY: it
#                              just rescales the axes in psmc_plot.pl, so every
#                              mu reuses the same .psmc model and only re-plots.
#
#   Net effect: N_patterns model runs + (N_patterns * N_mu) plots, with the
#   slow Steps 1-4 done a single time.
#
# DIRECTORY LAYOUT (unchanged; parameter-aware):
#   ${OUTDIR}/
#   ├── shared/                      built ONCE; independent of MU/GEN/pattern
#   │   ├── chrom_scaffolds.regions
#   │   ├── depth_cutoffs.txt
#   │   ├── <SAMPLE>.diploid.fq.gz
#   │   └── <SAMPLE>.psmcfa
#   ├── models/<PATTERN_TAG>/        one dir per PSMC_PATTERN (the model)
#   │   ├── <SAMPLE>.psmc
#   │   ├── <SAMPLE>.split.psmcfa
#   │   ├── boot_*.psmc
#   │   ├── <SAMPLE>.combined.psmc
#   │   └── model_params.txt
#   └── plots/<PATTERN_TAG>/mu<MU>_g<GEN>/   one dir per (pattern, MU) scaling
#       ├── <SAMPLE>_psmc_<PATTERN_TAG>_mu<MU>.{eps,pdf}
#       ├── <SAMPLE>_psmc_bootstrap_<PATTERN_TAG>_mu<MU>.{eps,pdf}
#       └── plot_params.txt
#
# One-time setup:
#   cd $my_softwares && git clone https://github.com/lh3/psmc
#   cd psmc && make && cd utils && make
#   (and: bash install_plotting_deps.sh   for gnuplot + ghostscript)
# =============================================================================

# ── Environment (matches your GONE pipeline) ─────────────────────────────────
source PG_project_env.sh
source ~/.bashrc

echo "Project Directory: $proj_dir"
echo "Data Directory   : $PG_data"
echo "Results Directory: $PG_results"

export PATH=${my_softwares}/envs/plotting/bin:${my_softwares}/bin:${my_softwares}/psmc:${my_softwares}/psmc/utils:$PATH
export LD_LIBRARY_PATH=${my_softwares}/lib:$LD_LIBRARY_PATH

ml samtools
ml bcftools

# ── USER SETTINGS — CONFIRM THESE ────────────────────────────────────────────
# 1. The female you chose (highest coverage + LOWEST heterozygosity/ROH)
SAMPLE="Hesperapis_oraria_2"
BAM="/90daydata/beenome100/hesperapis_oraria_genomics/Gulf-Coast-bee2/results/Long_read_Variant_Calling/Alignment/${SAMPLE}.bam"
REF="$ref_genome"   # <-- from PG_project_env.sh; confirm it's exported

# 2. Mutation rates to sweep (per bp per generation). PLOT-ONLY: each value
#    just rescales the axes. Add/remove freely; inner loop, no model rerun.
MU_LIST=("3.4e-9" "5.3e-9" "7.2e-9")

# 3. Generation time (years). PLOT-ONLY (univoltine -> g = 1).
#    Kept scalar; to also sweep g, make it GEN_LIST and add a loop under MU.
GEN="1"

# 4. Coverage for depth cutoffs (AUTO = measure from BAM).
MEANDEPTH="18"

# 5. PSMC atomic time-interval patterns to sweep. Each defines a MODEL dir;
#    outer loop, Steps 5-6 run once per pattern.
PSMC_PATTERN_LIST=("4+25*2+4+6" "1+1+1+1+25*2+4+6" "2+2+25*2+4+6")

# ── STEP TOGGLES ─────────────────────────────────────────────────────────────
# Steps 1-4 are global (run once). Steps 5-7 run inside the loops.
# Each step checks its OWN inputs, so skipping is safe as long as a later
# ENABLED step finds what it needs on disk. Common recipes:
#   New MU only ............. STEP1..6=false, STEP7=true   (reuses every model)
#   New pattern(s) only ..... STEP1..4=false, STEP5..7=true (reuses shared/)
#   Re-bootstrap+plot ....... STEP1..4=false, STEP5=false, STEP6,7=true
#   First full run .......... all true
RUN_STEP1=false   # restrict to 17 chromosome scaffolds  -> shared/        (once)
RUN_STEP2=false   # depth cutoffs                         -> shared/        (once)
RUN_STEP3=false   # diploid consensus FASTQ (SLOW)        -> shared/        (once)
RUN_STEP4=false   # fq2psmcfa (+ het check)               -> shared/        (once)
RUN_STEP5=false    # PSMC model run            -> models/<tag>/    (per pattern)
RUN_STEP6=false    # bootstrap                 -> models/<tag>/    (per pattern)
RUN_STEP7=true    # plot (uses MU, GEN)       -> plots/<tag>/mu_g/ (per pattern x mu)

# ── Analysis parameters ──────────────────────────────────────────────────────
OUTDIR="${PG_results}/PSMC/${SAMPLE}"
CHROM_SCAFFOLDS="scaffold0020 scaffold0002 scaffold0003 scaffold0004 scaffold0022 \
scaffold0006 scaffold0007 scaffold0023 scaffold0024 scaffold0032 scaffold0011 \
scaffold0043 scaffold0013 scaffold0045 scaffold0015 scaffold0016 scaffold0047"
MIN_MAPQ=30
MIN_BASEQ=20

PSMC_OPTS="-N25 -t15 -r5" # kept fixed; recorded in model_params.txt
NBOOT=100
THREADS=8

# ── Shared (pattern/MU-independent) paths ────────────────────────────────────
SHARED_DIR="${OUTDIR}/shared"
mkdir -p "${SHARED_DIR}"
REGIONS="${SHARED_DIR}/chrom_scaffolds.regions"
CUTOFFS="${SHARED_DIR}/depth_cutoffs.txt"
FQ="${SHARED_DIR}/${SAMPLE}.diploid.fq.gz"
PSMCFA="${SHARED_DIR}/${SAMPLE}.psmcfa"

echo "======================================================"
echo " PSMC : ${SAMPLE}"
echo " Patterns : ${PSMC_PATTERN_LIST[*]}"
echo " MU sweep : ${MU_LIST[*]}   (g=${GEN})"
echo " shared   : ${SHARED_DIR}"
echo " Started  : $(date)"
echo "======================================================"

# ── Step 0: sanity checks (always run) ───────────────────────────────────────
for f in "${BAM}" "${REF}"; do
    [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
done
[[ -s "${REF}.fai" ]] || samtools faidx "${REF}"
[[ -s "${BAM}.bai" ]] || samtools index "${BAM}"
command -v psmc        >/dev/null || { echo "ERROR: psmc not on PATH (build it, see header)"; exit 1; }
command -v vcfutils.pl >/dev/null || { echo "ERROR: vcfutils.pl not found (ships with bcftools misc)"; exit 1; }

# =============================================================================
# SHARED STEPS 1-4  (run once, before any loop)
# =============================================================================

# ── Step 1: restrict to the 17 chromosome-level scaffolds  -> shared/ ─────────
if [[ "${RUN_STEP1}" == true ]]; then
    > "${REGIONS}"
    for s in ${CHROM_SCAFFOLDS}; do
        len=$(awk -v s="$s" '$1==s {print $2}' "${REF}.fai")
        [[ -n "$len" ]] || { echo "ERROR: scaffold '$s' not in ${REF}.fai (check naming)"; exit 1; }
        printf "%s\t1\t%s\n" "$s" "$len" >> "${REGIONS}"
    done
    echo "[1] Chromosome scaffolds included: $(wc -l < "${REGIONS}")"
else
    echo "[1] SKIPPED (reuse ${REGIONS} if present)"
fi

# ── Step 2: depth cutoffs  -> shared/ ────────────────────────────────────────
if [[ "${RUN_STEP2}" == true ]]; then
    if [[ "${MEANDEPTH}" == "AUTO" ]]; then
        echo "[2] Estimating genome-wide mean depth (samtools coverage)..."
        MEANDEPTH=$(samtools coverage "${BAM}" \
            | awk 'NR>1 {len=$3-$2; sumd+=$7*len; sumlen+=len} END{printf "%.1f", sumd/sumlen}')
    fi
    MINDP=$(awk -v m="${MEANDEPTH}" 'BEGIN{d=int(m/3+0.5); if(d<1)d=1; print d}')
    MAXDP=$(awk -v m="${MEANDEPTH}" 'BEGIN{print int(2*m+0.5)}')
    echo "${MEANDEPTH} ${MINDP} ${MAXDP}" > "${CUTOFFS}"
    echo "[2] Mean depth ~ ${MEANDEPTH}x  ->  -d ${MINDP}  -D ${MAXDP}"
    (( $(awk -v m="${MEANDEPTH}" 'BEGIN{print (m<10)}') )) \
        && echo "    WARNING: mean depth < 10x; recent end of curve is unreliable."
elif [[ -s "${CUTOFFS}" ]]; then
    read MEANDEPTH MINDP MAXDP < "${CUTOFFS}"
    echo "[2] SKIPPED (reusing -d ${MINDP} -D ${MAXDP} from ${CUTOFFS})"
else
    echo "[2] SKIPPED (no cached cutoffs; fine unless Step 3 runs)"
fi

# ── Step 3: diploid consensus  -> shared/ ────────────────────────────────────
if [[ "${RUN_STEP3}" == true ]]; then
    [[ -n "${MINDP}" && -n "${MAXDP}" ]] \
        || { echo "ERROR: Step 3 needs depth cutoffs. Set RUN_STEP2=true once."; exit 1; }
    REGIONS_OPT=""
    if [[ -s "${REGIONS}" ]]; then REGIONS_OPT="-R ${REGIONS}"
    else echo "    WARNING: no ${REGIONS} -> mpileup uses ALL scaffolds. Run Step 1."; fi
    echo "[3] Building diploid consensus FASTQ -> ${FQ}"
    bcftools mpileup -q ${MIN_MAPQ} -Q ${MIN_BASEQ} \
            --max-depth $(( MAXDP*2 > 250 ? MAXDP*2 : 250 )) ${REGIONS_OPT} \
            -f "${REF}" "${BAM}" \
      | bcftools call -c \
      | vcfutils.pl vcf2fq -d ${MINDP} -D ${MAXDP} -Q ${MIN_BASEQ} \
      | gzip > "${FQ}"
else
    echo "[3] SKIPPED (reuse ${FQ})"
fi

# ── Step 4: fq2psmcfa  -> shared/ ────────────────────────────────────────────
if [[ "${RUN_STEP4}" == true ]]; then
    [[ -s "${FQ}" ]] || { echo "ERROR: ${FQ} not found. Enable Step 3."; exit 1; }
    echo "[4] fq2psmcfa -> ${PSMCFA}"
    fq2psmcfa -q20 "${FQ}" > "${PSMCFA}"
    HET=$(awk '!/^>/{for(i=1;i<=length($0);i++){c=substr($0,i,1); if(c=="K")h++; if(c!="N")t++}}
               END{if(t>0) printf "%.5f", h/t; else print "NA"}' "${PSMCFA}")
    echo "    Approx binned heterozygosity: ${HET}  (very low -> check coverage / inbreeding)"
else
    echo "[4] SKIPPED (reuse ${PSMCFA})"
fi

# =============================================================================
# OUTER LOOP over PSMC_PATTERN  (Steps 5-6 build the model per pattern)
#   INNER LOOP over MU           (Step 7 re-plots per mu, no model rerun)
# =============================================================================
for PSMC_PATTERN in "${PSMC_PATTERN_LIST[@]}"; do

    # PATTERN_TAG: make PSMC_PATTERN directory-safe ( * -> x , + -> _ )
    PATTERN_TAG="p$(printf '%s' "${PSMC_PATTERN}" | sed 's/\*/x/g; s/+/_/g')"
    MODEL_DIR="${OUTDIR}/models/${PATTERN_TAG}"
    PSMC_OUT="${MODEL_DIR}/${SAMPLE}.psmc"
    SPLIT="${MODEL_DIR}/${SAMPLE}.split.psmcfa"
    COMBINED="${MODEL_DIR}/${SAMPLE}.combined.psmc"
    mkdir -p "${MODEL_DIR}"

    echo ""
    echo "######################################################"
    echo "# PATTERN ${PSMC_PATTERN}   (tag ${PATTERN_TAG})"
    echo "# model : ${MODEL_DIR}"
    echo "######################################################"

    # ── Step 5: PSMC model run  -> models/<tag>/ ─────────────────────────────
    if [[ "${RUN_STEP5}" == true ]]; then
        [[ -s "${PSMCFA}" ]] || { echo "ERROR: ${PSMCFA} not found. Enable Step 4."; exit 1; }
        echo "[5][${PATTERN_TAG}] Running PSMC (pattern ${PSMC_PATTERN}) -> ${MODEL_DIR}"
        psmc ${PSMC_OPTS} -p "${PSMC_PATTERN}" -o "${PSMC_OUT}" "${PSMCFA}"
        { echo "PSMC_PATTERN=${PSMC_PATTERN}"
          echo "psmc_opts=${PSMC_OPTS}"
          echo "psmcfa=${PSMCFA}"
          echo "date=$(date)"; } > "${MODEL_DIR}/model_params.txt"
    else
        echo "[5][${PATTERN_TAG}] SKIPPED (reuse ${PSMC_OUT})"
    fi

    # ── Step 6: bootstrap  -> models/<tag>/ ──────────────────────────────────
    if [[ "${RUN_STEP6}" == true ]]; then
        [[ -s "${PSMCFA}" ]]   || { echo "ERROR: ${PSMCFA} not found. Enable Step 4."; exit 1; }
        [[ -s "${PSMC_OUT}" ]] || { echo "ERROR: ${PSMC_OUT} not found. Enable Step 5."; exit 1; }
        echo "[6][${PATTERN_TAG}] Bootstrapping ${NBOOT} reps (pattern ${PSMC_PATTERN})..."
        splitfa "${PSMCFA}" > "${SPLIT}"
        seq ${NBOOT} | xargs -P ${THREADS} -I{} \
            psmc ${PSMC_OPTS} -b -p "${PSMC_PATTERN}" \
                 -o "${MODEL_DIR}/boot_{}.psmc" "${SPLIT}"
        cat "${PSMC_OUT}" "${MODEL_DIR}"/boot_*.psmc > "${COMBINED}"
    else
        echo "[6][${PATTERN_TAG}] SKIPPED (reuse ${COMBINED})"
    fi

    # ── Step 7: plot per MU  -> plots/<tag>/mu<MU>_g<GEN>/ ────────────────────
    # psmc_plot.pl draws an EPS via gnuplot; convert EPS->PDF with ps2pdf (so the
    # -p/epstopdf path is intentionally NOT used). Tools from install_plotting_deps.sh.
    if [[ "${RUN_STEP7}" == true ]]; then
        [[ -s "${PSMC_OUT}" ]] \
            || { echo "ERROR: ${PSMC_OUT} not found. Run Step 5 for pattern ${PSMC_PATTERN} first."; exit 1; }
        command -v gnuplot >/dev/null \
            || { echo "ERROR: gnuplot not on PATH. Run install_plotting_deps.sh first."; exit 1; }

        for MU in "${MU_LIST[@]}"; do
            PLOT_DIR="${OUTDIR}/plots/${PATTERN_TAG}/mu${MU}_g${GEN}"
            mkdir -p "${PLOT_DIR}"
            echo "[7][${PATTERN_TAG}] Plotting (u=${MU}, g=${GEN}) -> ${PLOT_DIR}"

            # Run plotting in a subshell so the cd does not leak across iterations.
            (
                cd "${PLOT_DIR}"
                { echo "MU=${MU}"; echo "GEN=${GEN}"; echo "PSMC_PATTERN=${PSMC_PATTERN}"
                  echo "model_dir=${MODEL_DIR}"; echo "date=$(date)"; } > plot_params.txt

                plot_psmc() {   # $1 = output prefix (relative) , $2 = input .psmc (absolute)
                    psmc_plot.pl -u "${MU}" -g "${GEN}" "$1" "$2"
                    if command -v ps2pdf >/dev/null && [[ -s "$1.eps" ]]; then
                        ps2pdf -dEPSCrop "$1.eps" "$1.pdf" && echo "    wrote ${PLOT_DIR}/$1.pdf"
                    else
                        echo "    wrote ${PLOT_DIR}/$1.eps (ps2pdf unavailable -> EPS only)"
                    fi
                }

                plot_psmc "${SAMPLE}_psmc_${PATTERN_TAG}_mu${MU}" "${PSMC_OUT}"
                if [[ -s "${COMBINED}" ]]; then
                    plot_psmc "${SAMPLE}_psmc_bootstrap_${PATTERN_TAG}_mu${MU}" "${COMBINED}"
                else
                    echo "    (no combined.psmc for pattern ${PSMC_PATTERN} -> skip bootstrap plot; enable Step 6)"
                fi
            )
        done
    else
        echo "[7][${PATTERN_TAG}] SKIPPED (plotting)"
    fi

done

echo ""
echo "======================================================"
echo " PSMC complete : $(date)"
echo " shared inputs : ${SHARED_DIR}/"
echo " models        : ${OUTDIR}/models/    (one dir per pattern)"
echo " plots         : ${OUTDIR}/plots/     (one dir per pattern x mu)"
echo "======================================================"


# =============================================================================
# EXPORT: gather the combined .psmc models for off-cluster transfer
#   .psmc is mutation-rate INDEPENDENT (mu is applied only at plotting time),
#   so ONE combined.psmc per PATTERN regenerates every (pattern x mu) plot.
#   Files are tagged by PATTERN_TAG so the export dir enumerates all combos.
# =============================================================================
RUN_EXPORT=true
if [[ "${RUN_EXPORT}" == true ]]; then
    EXPORT_DIR="${OUTDIR}/plots/export_combined_psmc"
    mkdir -p "${EXPORT_DIR}"
    echo ""
    echo "[export] Collecting combined .psmc models -> ${EXPORT_DIR}"

    n_copied=0
    for PSMC_PATTERN in "${PSMC_PATTERN_LIST[@]}"; do
        PATTERN_TAG="p$(printf '%s' "${PSMC_PATTERN}" | sed 's/\*/x/g; s/+/_/g')"
        SRC="${OUTDIR}/models/${PATTERN_TAG}/${SAMPLE}.combined.psmc"
        PARAMS="${OUTDIR}/models/${PATTERN_TAG}/model_params.txt"
        if [[ -s "${SRC}" ]]; then
            cp -p "${SRC}" "${EXPORT_DIR}/${SAMPLE}.${PATTERN_TAG}.combined.psmc"
            [[ -s "${PARAMS}" ]] && \
                cp -p "${PARAMS}" "${EXPORT_DIR}/${SAMPLE}.${PATTERN_TAG}.model_params.txt"
            runs=$(grep -c '^//' "${SRC}" 2>/dev/null || echo '?')
            echo "    + ${PATTERN_TAG}  (${runs} runs; expect 101 = 1 main + 100 boot)"
            ((n_copied++))
        else
            echo "    ! ${PATTERN_TAG}: ${SRC} missing (run Steps 5-6) -> skipped"
        fi
    done

    if (( n_copied > 0 )); then
        TARBALL="${OUTDIR}/plots/${SAMPLE}_combined_psmc_export.tar.gz"
        tar -czf "${TARBALL}" -C "${EXPORT_DIR}" .
        echo "[export] ${n_copied} model(s) collected"
        echo "[export] tarball : ${TARBALL}"
        echo "[export] fetch e.g.: scp <user>@ceres.scinet.usda.gov:${TARBALL} ."
    else
        echo "[export] nothing copied — no combined.psmc found. Did Step 6 run?"
    fi
fi