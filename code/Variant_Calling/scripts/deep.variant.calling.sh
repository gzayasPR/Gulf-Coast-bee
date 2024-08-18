#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --job-name="$1"
#SBATCH --output=$1.out
#SBATCH --error=$1.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=72:00:00  # Adjusted walltime to 24 hours
#SBATCH --ntasks=48  # Using all available cores on a single node
#SBATCH --nodes=1  # Use a single node
#SBATCH --partition=atlas  # Selecting the atlas partition
#SBATCH --qos=normal  # Using the normal QoS


name=$1
bam_dir=$2
ref_genome=$3
known_vcf=$4
ploidy=$5
output_dir=$6/DeepVariant/
deep_variant=/project/90daydata/beenome100/hesperapis_oraria_genomics/Hesperapis_oraria_Pop_genomics/code/DeepVariant
mkdir -p ${output_dir}
ml miniconda3
ml apptainer
eval "$(conda shell.bash hook)"
input_bam=${bam_dir}/${name}.marked_duplicates.bam

output_vcf=${output_dir}/${name}.deepvariant.vcf.gz
output_gvcf=${output_dir}/${name}.deepvariant.g.vcf.gz
cd $deep_variant
export APPTAINER_CACHEDIR=$TMPDIR 
export APPTAINER_TMPDIR=$TMPDIR
scaffolds=$(awk '{printf "%s%s",sep,$1; sep=","} END{print ""}' ${ref_genome}.fai)
if [ "$ploidy" == "1" ];
then
echo "PLoidy 1"
bed_file=${output_dir}/${name}.bed
awk '{print $1"\t0\t"$2}' ${ref_genome}.fai > ${bed_file}
    # Run DeepVariant using Apptainer
    apptainer exec \
      --bind ${bam_dir}:/input \
      --bind ${output_dir}:/output \
      ${deep_variant}/deepvariant_1.6.1.sif \
      /opt/deepvariant/bin/run_deepvariant \
      --model_type=WGS \
      --haploid_contigs=$scaffolds \
      --ref=${ref_genome} \
      --reads=${input_bam} \
      --output_vcf=${output_vcf} \
      --output_gvcf=${output_gvcf}  \ 
     --num_shards=48 \
     --par_regions_bed ${bed_file} \
      --intermediate_results_dir "${output_dir}/intermediate_results_dir" 
       fi
if [ "$ploidy" == "2" ];
then
echo "PLoidy 2"
apptainer exec \
  --bind ${bam_dir}:/input \
  --bind ${output_dir}:/output \
  ${deep_variant}/deepvariant_1.6.1.sif \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=${ref_genome} \
  --reads=${input_bam} \
  --output_vcf=${output_vcf} \
  --output_gvcf=${output_gvcf}  \ 
     --num_shards=48 \
  --intermediate_results_dir "${output_dir}/intermediate_results_dir" 
  fi