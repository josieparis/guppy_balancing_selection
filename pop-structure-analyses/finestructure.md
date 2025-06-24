### These are series of SLURM scripts
### Except script 3 where stuff gets complicated ... 

Script 1: 01_convertvcf2haps.sh

```
#!/bin/bash
#SBATCH -D .
#SBATCH -p sq
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=vcf2haps
#SBATCH --error=logs/vcf2haps.err.txt
#SBATCH --output=logs/vcf2haps.out.txt
#SBATCH --export=All

# This assumes that all separate chromosomes have been phased and are in a directory somewhere...

module purge
module load PLINK/2.00a2.3_x86_64 

PHASED_VCF=<dir>/holi11_vcf_files/phased_vcfs
DATASET=holi11
OUT=<dir>/NFDS_analysis/pop_structure/finestructure/data_holi11

## loop over chromosomes and export haplotypes
for chr in $(cut -f1 ~/startup/STAR/STAR.chromosomes.release.fasta.fai | grep -v "alt" | grep -v "chr12")
do
plink2 --vcf ${PHASED_VCF}/${DATASET}_phased_${chr}_shapeit_phased.gt.vcf.gz --export haps --out ${OUT}/${DATASET}_phased_${chr}_shapeit_phased \
--allow-extra-chr --double-id
done
```

Script 2: 02_finestructure_convert_shapeit2chromo.sh

```
#!/bin/bash
#SBATCH -D .
#SBATCH -p highmem
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=conversion_for_fineSTRUC
#SBATCH --error=logs/conversion_for_fineSTRUCTURE.err.txt
#SBATCH --output=logs/conversion_for_fineSTRUCTURE.out.txt
#SBATCH --export=All

# This assumes that all separate chromosomes have been phased with shapeit and are in a directory somewhere...

module load Perl/5.28.0-GCCcore-7.3.0

MASTER=/lustre/home/jrp228/NERC/people/josie/NFDS_analysis/pop_structure/finestructure
DATASET=holi11

# Phased chromosomes can be cat together to give full dataset
for chr in $(cut -f1 ~/startup/STAR/STAR.chromosomes.release.fasta.fai | grep -v "alt" | grep -v "chr12")
do
cat $MASTER/data_holi11/*_${chr}_*.haps >> $MASTER/data_holi11/${DATASET}_allchr_shapeit_beagle_phase.haps
done

# Now convert to chromopainter
perl ~/programs/finestructure_4.1.1/impute2chromopainter.pl -J ${MASTER}/data_holi11/${DATASET}_allchr_shapeit_beagle_phase.haps ${MASTER}/data_holi11/${DATASET}_chromoIN

# Make rec file
perl ~/programs/finestructure_4.1.1/makeuniformrecfile.pl ${MASTER}/data_holi11/${DATASET}_chromoIN.phase ${MASTER}/data_holi11/${DATASET}_chromoIN.rec
```

Script 3: 03_run_fineSTRUCTURE_jrp.sh

```
#!/bin/bash

# -----------------------
# NOTE - THIS SCRIPT IS RUN FROM THE COMMANDLINE, THE SUBMISSION OCCURS THROUGH AN EDITED VERSION OF FINESTRUCTURE'S qsub_run.sh script
# -----------------------

# Set up environment
MASTER=<dir>NERC/people/josie/NFDS_analysis/pop_structure/finestructure
FINE_S=<dir>/programs/finestructure_4.1.1
DATASET="holi11"

# I have added this so that the libgsl.so.0 symlink is in the library path. Links to ~./linuxbrew/libs/libgsl.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/home/jrp228/programs/finestructure_4.1.1

# Load modules
module load GLib/2.53.5-GCCcore-6.4.0
module load GSL/2.4-GCCcore-6.4.0
module load Perl/5.28.0-GCCcore-7.3.0
module load R/3.5.1-foss-2018b

echo '**** START PAINT/fineSTRUCTURE ****'

# remove duplicates from ids file
#uniq $MASTER/outputs/${DATASET}_chromoIN.ids > $MASTER/outputs/${DATASET}_chromoIN.uniq.ids

# Change
#cd $MASTER/data

# -------------------------------------
# STAGES 1-2 TAKE APPROX 12 HOURS ON PQ
# -------------------------------------

echo '**** Write Command Files 1 ****'

$FINE_S/fs $MASTER/outputs_holi11/${DATASET}_fineSTRUCTURE.cp \
-hpc 1 \
-idfile $MASTER/outputs_holi11/${DATASET}_chromoIN.ids \
-phasefiles $MASTER/outputs_holi11/${DATASET}_chromoIN.phase \
-recombfiles $MASTER/outputs_holi11/${DATASET}_chromoIN.rec \
-s1emits 10 -s1minsnps 10000 -go

# We need to edit the commandfiles to replace the fs function with where we store the new fs command
#sed -i 's/fs /fineSTRUCTURE /g' $MASTER/outputs_holi11/${DATASET}_fineSTRUCTURE/commandfiles/commandfile1.txt

# We have to edit the path to the outputs because it does weird things ...

echo '**** Run Command Files 1 on HPC ****'
# NOTES -
# I have edited the qsub script a fair bit.
# The script must be run from the commandfiles directory, this is so the qsub_run script includes all the correct paths for outputs
# However the qsub_run script now includes a cd ../../ which brings the actual analysis back to the outputs directory so that all of the inputs are available (these do not have full paths!)
cd $MASTER/outputs_holi11/${DATASET}_fineSTRUCTURE/commandfiles/


######################
# RUN STAGE 1 ON COMMAND LINE
######################
#$FINE_S/qsub_run.sh -f commandfile1.txt -n 16 -m 16 -w 3 -P parallel # 16 ppn per node, 16 jobs per node, 3 hours walltime
$FINE_S/sbatch_run.sh -q sbatch -f commandfile1.txt -n 4 -m 4 -w 3 -P parallel # 16 ppn per node, 16 jobs per node, 3 hours walltime
######################
echo '**** Finished Command Files 1 on HPC ****'

cd $MASTER/outputs_holi11

# Stage 1 successful, reset to stage 2 (Need to change the .cp file paths!)
# Need to run a number of processes that doesn't overwhelm the RAM limitations. Inds per process of 20 yields 6 jobs which keeps RAM in check
$FINE_S/fs ${DATASET}_fineSTRUCTURE.cp -indsperproc 20 -go

# Reset if needs be
#fineSTRUCTURE ${DATASET}_fineSTRUCTURE.cp -reset 2 -go

# Edit the command files again
# sed -i 's/fs /fineSTRUCTURE /g'  $MASTER/outputs_holi11/${DATASET}_fineSTRUCTURE/commandfiles/commandfile2.txt

cd ${DATASET}_fineSTRUCTURE/commandfiles/
######################
# RUN ON STAGE 2 COMMAND LINE - TOTAL OF 9 JOBS TO RUN (INDS PER PROC = 20) Ran in ~72 hours
######################
$FINE_S/sbatch_run_highmem.sh -f commandfile2.txt -n 16 -m 9 -w 80 -P parallel # 16 ppn per node, 9 jobs per node, 24 hours walltime
######################

cd $MASTER/outputs_holi11

$FINE_S/fs holi11_fineSTRUCTURE.cp -reset 2 -go
# holi11 - Successfully run ChromoCombine stage! Inferred a 'c' value of 0.233028

cd $MASTER/outputs_holi11

#Stage 2 successful, reset to stage 3 (runs very quickly)
$FINE_S/fs ${DATASET}_fineSTRUCTURE.cp -reset 3 -go
$FINE_S/fs  ${DATASET}_fineSTRUCTURE.cp -s3iters 1000000 -allowdep 1 -hpc 1 -go

# Edit and move
sed -i 's/fs/fineSTRUCTURE/g'  $MASTER/outputs_holi11/${DATASET}_fineSTRUCTURE/commandfiles/commandfile3.txt
cd ${DATASET}_fineSTRUCTURE/commandfiles/

######################
# RUN ON STAGE 3 COMMAND LINE - TOTAL OF 9 JOBS TO RUN (INDS PER PROC = 20) runs fast >6 hrs
######################
${MASTER}/scripts/qsub_run_highmem.sh -f commandfile3.txt -n 16 -w 96
######################

# holi_13 = Conversion looks good

# Stage 3 successful, reset to stage 4 (runs very quickly)
fineSTRUCTURE ${DATASET}_fineSTRUCTURE.cp -reset 4 -go

#Need to edit commandfile4 so that 'fs fs' becomes 'fineSTRUCTURE fs'
sed -i 's/fs/fineSTRUCTURE/' $MASTER/outputs/${DATASET}_fineSTRUCTURE/commandfiles/commandfile4.txt

cd ${DATASET}_fineSTRUCTURE/commandfiles/
######################
# RUN ON STAGE 4 COMMAND LINE - I sometimes just run this in an interactive job, it is quick. NOTE - run from outputs/ obviously...
######################
${MASTER}/scripts/qsub_run.sh -f commandfile4.txt -n 16 -m 1 -w 6 -p
######################
echo '**** END PAINT/fineSTRUCTURE ****'

######################
# PLOT OUTPUT USING R
######################
Rscript ${MASTER}/R/Plot_fineSTRUCTURE_tree.R ${DATASET}
```









