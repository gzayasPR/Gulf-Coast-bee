#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=Fst_ANGSD_%j.out
#SBATCH --error=Fst_ANGSD_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb

source ANGSD_project_env.sh
source ~/.bashrc
# Set the correct environment variables
export PATH=${my_softwares}/bin:$PATH
export LD_LIBRARY_PATH=${my_softwares}/lib:$LD_LIBRARY_PATH
export CFLAGS="-I${my_softwares}/include"
export LDFLAGS="-L${my_softwares}/lib"
# Define project directories
echo "Project Directory: $proj_dir"
echo "Code Directory: $ANGSD_code"
echo "Data Directory: $ANGSD_data"
echo "Results Directory: $ANGSD_results"


out_dir=${ANGSD_results}/Fst/
mkdir -p $out_dir

angsd_dir=${my_softwares}/angsd/env
angsd_beagle=${ANGSD_results}/ANGSD_out/combined/combined.beagle.gz
angsd_beagle_mafs=${ANGSD_results}/ANGSD_out/combined/combined.mafs.gz
pop_info=${ANGSD_results}/ANGSD_out/combined.pop.info
bam_dir=${proj_dir}/results/Variant_Calling/3.Cleaning/
bam_filelist=${ANGSD_results}/ANGSD_out/combined_bam.filelist

cd ${out_dir}
zcat ${angsd_beagle_mafs} | cut -f 1,2 | tail -n +2 > ${out_dir}/sites.txt

source ~/.bashrc
ml miniconda3
conda activate $angsd_dir
angsd sites index ${out_dir}/sites.txt
out_dir_pop2=${ANGSD_results}/Fst/ALxFL/
ml samtools
# Index the reference genome if not already indexed
samtools faidx ${ref_genome}
mkdir -p $out_dir_pop2
grep "USA-FL_Escambia_co" ${meta_data} | awk -F "," '{print $1}' > ${out_dir_pop2}/pop1.names

grep -f ${out_dir_pop2}/pop1.names ${bam_filelist} > ${out_dir_pop2}/pop1.bam.filelist

grep "USA-AL_Baldwin_co" ${meta_data} | awk -F "," '{print $1}' > ${out_dir_pop2}/pop2.names
grep -f ${out_dir_pop2}/pop2.names ${bam_filelist} > ${out_dir_pop2}/pop2.bam.filelist


#this is with 2pops
#first calculate per pop saf for each populatoin
angsd -b ${out_dir_pop2}/pop1.bam.filelist  -GL 1 -doSaf 1  -anc ${ref_genome} -sites ${out_dir}/sites.txt  -doCounts 1 -minQ 20 -minMapQ 30 -remove_bads 1 -only_proper_pairs 1 -P 8 -out  ${out_dir_pop2}/pop1 
angsd -b ${out_dir_pop2}/pop2.bam.filelist  -GL 1 -doSaf 1 -anc ${ref_genome} -sites ${out_dir}/sites.txt  -doCounts 1 -minQ 20 -minMapQ 30 -remove_bads 1 -only_proper_pairs 1   -P 8 -out  ${out_dir_pop2}/pop2 
#calculate the 2dsfs prior
realSFS -P 8 ${out_dir_pop2}/pop1.saf.idx ${out_dir_pop2}/pop2.saf.idx > ${out_dir_pop2}/pop1.pop2.2dsfs
#prepare the fst for easy window analysis etc
realSFS fst index ${out_dir_pop2}/pop1.saf.idx ${out_dir_pop2}/pop2.saf.idx -sfs ${out_dir_pop2}/pop1.pop2.2dsfs -fstout ${out_dir_pop2}/here -P 8
#get the global estimate
realSFS  fst stats ${out_dir_pop2}/here.fst.idx -P 8 > ${out_dir_pop2}/pop1_pop2_fst.stats
#below is not tested that much, but seems to work
realSFS  fst stats2 ${out_dir_pop2}/here.fst.idx -win 50000 -step 10000 > ${out_dir_pop2}/slidingwindow


cd ${out_dir}
out_dir_pop3=${ANGSD_results}/Fst/admix_pops/
mkdir -p $out_dir_pop3
meta_data=${ANGSD_results}/ADMIXTURE/admix.csv
grep "Q1" ${meta_data} | awk -F "," '{print $1}' > ${out_dir_pop3}/pop1.names
grep -f ${out_dir_pop3}/pop1.names ${bam_filelist} > ${out_dir_pop3}/pop1.bam.filelist
echo "Pop1"
cat ${out_dir_pop3}/pop1.bam.filelist

grep "Q2" ${meta_data} | awk -F "," '{print $1}' > ${out_dir_pop3}/pop2.names
grep -f ${out_dir_pop3}/pop2.names ${bam_filelist} > ${out_dir_pop3}/pop2.bam.filelist
echo "Pop2"
cat ${out_dir_pop3}/pop2.bam.filelist

grep "Q3" ${meta_data} | awk -F "," '{print $1}' > ${out_dir_pop3}/pop3.names
grep -f ${out_dir_pop3}/pop3.names ${bam_filelist} > ${out_dir_pop3}/pop3.bam.filelist
echo "Pop3"
cat ${out_dir_pop3}/pop3.bam.filelist
#this is with 2pops
#first calculate per pop saf for each populatoin
angsd -b ${out_dir_pop3}/pop1.bam.filelist  -GL 1 -doSaf 1  -anc ${ref_genome} -sites ${out_dir}/sites.txt  -doCounts 1 -minQ 20 -minMapQ 30 -remove_bads 1 -only_proper_pairs 1 -P 8 -out  ${out_dir_pop3}/pop1 
angsd -b ${out_dir_pop3}/pop2.bam.filelist  -GL 1 -doSaf 1 -anc ${ref_genome} -sites ${out_dir}/sites.txt  -doCounts 1 -minQ 20 -minMapQ 30 -remove_bads 1 -only_proper_pairs 1   -P 8 -out  ${out_dir_pop3}/pop2
angsd -b ${out_dir_pop3}/pop3.bam.filelist  -GL 1 -doSaf 1 -anc ${ref_genome} -sites ${out_dir}/sites.txt  -doCounts 1 -minQ 20 -minMapQ 30 -remove_bads 1 -only_proper_pairs 1   -P 8 -out  ${out_dir_pop3}/pop3
#calculate all pairwise 2dsfs's
realSFS -P 8 ${out_dir_pop3}/pop1.saf.idx ${out_dir_pop3}/pop2.saf.idx > ${out_dir_pop3}/pop1.pop2.2dsfs
realSFS -P 8 ${out_dir_pop3}/pop1.saf.idx ${out_dir_pop3}/pop3.saf.idx > ${out_dir_pop3}/pop1.pop3.2dsfs
realSFS -P 8 ${out_dir_pop3}/pop2.saf.idx ${out_dir_pop3}/pop3.saf.idx > ${out_dir_pop3}/pop2.pop3.2dsfs

realSFS fst index ${out_dir_pop3}/pop1.saf.idx ${out_dir_pop3}/pop2.saf.idx ${out_dir_pop3}/pop3.saf.idx \
        -sfs ${out_dir_pop3}/pop1.pop2.2dsfs \
        -sfs ${out_dir_pop3}/pop1.pop3.2dsfs \
        -sfs ${out_dir_pop3}/pop2.pop3.2dsfs \
        -fstout ${out_dir_pop3}/here -P 8
#get the global estimate
realSFS  fst stats ${out_dir_pop3}/here.fst.idx -P 8 > ${out_dir_pop3}/pop1_pop2_pop3_fst.stats
#below is not tested that much, but seems to work
realSFS  fst stats2 ${out_dir_pop3}/here.fst.idx -win 50000 -step 10000 > ${out_dir_pop3}/slidingwindow

