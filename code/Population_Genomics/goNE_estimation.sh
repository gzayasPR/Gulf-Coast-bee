#!/bin/bash
#SBATCH --job-name=gone_ne
#SBATCH --output=gone_ne_%j.out
#SBATCH --error=gone_ne_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your@email.edu

# =============================================================================
# GONE Ne Estimation Pipeline
#
# Estimates effective population size (Ne) over recent generations using
# linkage disequilibrium. Uses GONE (Santiago et al. 2020, MBE).
#
# Input  : PLINK .bed/.bim/.fam (your existing bee_hesperapis_plink files)
# Output : OUTPUT_Ne_* files in the run directory
#
# GONE GitHub : https://github.com/esrud/GONE
# Citation    : Santiago et al. 2020, Mol Biol Evol 37(12):3642-3653
#
# Setup (run once before submitting):
#   git clone https://github.com/esrud/GONE.git
#   chmod +x GONE/Linux/GONE
#   chmod +x GONE/Linux/script_GONE.sh
# =============================================================================

# ── User settings ─────────────────────────────────────────────────────────────



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
out_dir=${PG_results}/GONE/${variant_caller}
vcf_file=${PG_data}/variants/${variant_caller}/Hesperapis_oraria.vcf
meta_data=${PG_data}/variants/Samples.Metadata.csv
name=Hesperapis_oraria_Hetero
king=${my_softwares}/king
ngsRelate_dir=${my_softwares}/ngsRelate
mkdir -p ${out_dir}
cd ${out_dir}
ml vcftools
ml bcftools 

min_Q=20
min_meanDP=3 # Setting the minimum read depth to 10 as per the new criteria
max_meanDP=55
max_locus_missing="0.95"
GQ_threshold=20
minDP=10 #Setting the minimum read depth to 10 as per the new criteria
maxDP=100

mkdir -p ${out_dir}
ml vcftools
ml bcftools
bcftools +fixploidy ${vcf_file}  -- -f 2 > ${out_dir}/diploid.vcf
# Create site-specific VCF files

# Filter each VCF file
vcftools --vcf ${out_dir}/diploid.vcf \
    --min-meanDP ${min_meanDP} \
    --max-meanDP ${max_meanDP} \
    --minQ ${min_Q} \
    --max-missing ${max_locus_missing} \
    --min-alleles 2 \
    --max-alleles 2 \
    --remove-indels \
    --out ${out_dir}/filtered \
    --recode --recode-INFO-all
            
ml plink2
plink2 -vcf ${out_dir}/filtered.recode.vcf -make-bed \
        --vcf-half-call h -out  ${out_dir}/bee_hesperapis_plink



PLINK_PREFIX="${out_dir}/bee_hesperapis_plink"   # no extension
GONE_DIR="$my_softwares/GONE/Linux/"                       # directory with GONE binary + script_GONE.sh
WORK_DIR="$out_dir/GONE_results"                            # output directory
SPECIES_NAME="Hesperapis_oraria"


# GONE INPUT_PARAMETERS
# PHASE      : 0 = unphased diploid, 1 = phased
# DIST       : 1 = Haldane (recommended for unphased data)
# NGEN       : max generations back to infer (200 = default)
# NBIN       : number of LD bins (400 = default)
# MINGAPSIZE : min gap between SNPs in bp
# MAXDIST    : max physical distance between SNP pairs in Mb
# HAPLOID    : 0 = diploid, 1 = haploid (male bees are haploid; set 1 if males included)
# MAF        : minimum minor allele frequency
# ZERO       : 0 = include zero LD values
# RUNID      : label for output files
write_params() {
cat > INPUT_PARAMETERS_FILE << 'PARAMS'
PHASE=0
cMMb=1
DIST=1
NGEN 2000
NBIN 400
MAF=0.00
ZERO=1
maxNCHROM=-99
maxNSHP=50000
hc=0.05
REPS=40
threads=-99
PARAMS
}
 
# Setup
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}" 
 
echo "======================================================"
echo " GONE Ne Estimation"
echo " Species : ${SPECIES_NAME}"
echo " PLINK   : ${PLINK_PREFIX}"
echo " Started : $(date)"
echo "======================================================"
 
# Step 1: Convert PLINK bed to ped/map
# GONE requires text ped/map format, not binary bed/bim/fam
# PLINK2 uses --export ped (not --recode as in PLINK1)
echo ""
echo "[1/4] Converting PLINK binary to ped/map..."
plink2 \
    --bfile "${PLINK_PREFIX}" \
    --allow-extra-chr \
    --export ped \
    --out "${SPECIES_NAME}"
 
# GONE requires integer chromosome IDs, not scaffold names
# Remap scaffold0001 etc. to sequential integers via the .map file
echo "[1b] Remapping scaffold names to integers for GONE..."
awk '{print $1}' "${SPECIES_NAME}.map" | sort -u > scaffold_list.txt
awk 'NR==FNR{map[$1]=NR; next} {$1=map[$1]; print}' \
    scaffold_list.txt "${SPECIES_NAME}.map" > "${SPECIES_NAME}_renum.map"
cp "${SPECIES_NAME}.ped" "${SPECIES_NAME}_renum.ped"
echo "  Scaffolds remapped: $(wc -l < scaffold_list.txt)"
 
# Step 2: Filter SNPs for GONE
# GONE works best with 10k-50k well-distributed SNPs
# Filter by MAF, missingness, and LD pruning
echo ""
echo "[2/4] Filtering SNPs (MAF > 0.05, LD pruning)..."
 
# Get prune-in list using PLINK2
plink2 \
    --ped "${SPECIES_NAME}_renum.ped" \
    --map "${SPECIES_NAME}_renum.map" \
    --allow-extra-chr \
    --maf 0.05 \
    --geno 0.1 \
    --bad-ld \
    --indep-pairwise 50 10 0.5 \
    --out "${SPECIES_NAME}_pruned"
 
# Extract pruned SNPs and export final ped/map for GONE
plink2 \
    --ped "${SPECIES_NAME}_renum.ped" \
    --map "${SPECIES_NAME}_renum.map" \
    --allow-extra-chr \
    --extract "${SPECIES_NAME}_pruned.prune.in" \
    --export ped \
    --out "${SPECIES_NAME}_gone_input"
 
echo "  SNPs retained: $(wc -l < ${SPECIES_NAME}_gone_input.map)"
echo "  Individuals  : $(wc -l < ${SPECIES_NAME}_gone_input.ped)"
 
# Step 3: Run GONE
echo ""
echo "[3/4] Running GONE..."
write_params
 
# GONE requires its binary and script_GONE.sh to be in the working directory
cp "${GONE_DIR}/INPUT_PARAMETERS_FILE" .
cp -r "${GONE_DIR}/PROGRAMMES" .
cp "${GONE_DIR}/script_GONE.sh" .
chmod a+rwx */*

# Run - argument is ped/map base name (no extension)
bash script_GONE.sh "${SPECIES_NAME}_gone_input"

# Step 4: Plot Ne trajectory in R
echo ""
echo "[4/4] Plotting Ne trajectory..."
Rscript - << 'RSCRIPT'
library(ggplot2)
 
# GONE output file: Output_Ne_<RUNID>
# Columns: Generation  Ne (no header)
ne_files <- list.files(".", pattern = "^Output_Ne_", full.names = TRUE)
if (length(ne_files) == 0) stop("No Output_Ne_ files found.")
 
ne <- read.table(ne_files[1], header = FALSE,
                 col.names = c("Generation", "Ne"))
ne <- ne[ne$Generation > 0 & ne$Ne > 0, ]
ne <- ne[order(ne$Generation), ]
 
# Report Ne at generations corresponding to ROH length thresholds
# Formula: ROH_Mb = 100 / (2 * g * cM_per_Mb)
# For bee (~20 cM/Mb): g = 100 / (2 * ROH_Mb * 20) = 2.5 / ROH_Mb
cat("\nNe at key generations:\n")
key_gens <- c(3, 5, 10, 20, 50, 100, 200)
ne_key   <- ne[ne$Generation %in% key_gens, ]
print(ne_key)
 
# Find sharpest Ne decline (bottleneck) as a natural threshold anchor
ne$delta <- c(NA, diff(ne$Ne))
bottleneck_gen <- ne$Generation[which.min(ne$delta)]
bottleneck_Ne  <- ne$Ne[ne$Generation == bottleneck_gen]
cat(sprintf("\nSharpest Ne decline at generation : %d (Ne = %.0f)\n",
            bottleneck_gen, bottleneck_Ne))
cat(sprintf("Equivalent ROH length at 1 cM/Mb  : %.2f Mb\n",
            100 / (2 * bottleneck_gen)))
cat(sprintf("Equivalent ROH length at 20 cM/Mb : %.3f Mb\n",
            100 / (2 * bottleneck_gen * 20)))
 
# Plot Ne trajectory
p <- ggplot(ne, aes(x = Generation, y = Ne)) +
  geom_line(colour = "#1A6B7C", linewidth = 0.9) +
  geom_point(size = 1.5, colour = "#1A6B7C") +
  geom_vline(xintercept = bottleneck_gen, linetype = "dashed",
             colour = "#E05C2A", linewidth = 0.7) +
  annotate("text",
           x     = bottleneck_gen * 1.1,
           y     = max(ne$Ne) * 0.9,
           label = paste0("Bottleneck\nGen ", bottleneck_gen),
           hjust = 0, size = 3.2, colour = "#E05C2A") +
  scale_x_log10(name = "Generations Ago") +
  scale_y_log10(name = "Effective Population Size (Ne)",
                labels = scales::comma) +
  labs(title    = "GONE: Ne Trajectory",
       subtitle = "Hesperapis oraria") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
 
ggsave("Ne_trajectory.pdf", p, width = 7, height = 5)
cat("Saved: Ne_trajectory.pdf\n")
 
# Save Ne table for use in ROH threshold decisions
write.table(ne[, c("Generation","Ne")], "Ne_summary.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")
cat("Saved: Ne_summary.txt\n")
RSCRIPT
 
echo ""
echo "======================================================"
echo " GONE complete : $(date)"
echo " Results in    : ${WORK_DIR}/"
echo ""
echo " Key output files:"
echo "   Output_Ne_*        Raw Ne per generation (from GONE)"
echo "   Ne_trajectory.pdf  Ne plot with bottleneck annotation"
echo "   Ne_summary.txt     Ne table for all generations"
echo ""
echo " Next step:"
echo "   Use Ne_summary.txt + roh_breakpoints.R output to set"
echo "   species-appropriate thresholds in plot_roh_ideogram_final.R"
echo "   (THR_INTERMEDIATE and THR_RECENT variables)"
echo "======================================================"