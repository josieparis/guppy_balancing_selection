#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16 # nodes=number of nodes required. ppn=number of processors per node
#SBATCH -A Research_Project-T109423
#SBATCH -p pq
#SBATCH --job-name=APHPballlermixchr
#SBATCH --error=APHPballermixchr.err.txt
#SBATCH --output=APHPballermixchr.out.txt
#SBATCH --export=All
#SBATCH -D .
#SBATCH --array=1-23%23

module load SciPy-bundle/2020.03-foss-2020a-Python-3.8.2
INPUT_DIR=/gpfs/ts0/projects/Research_Project-T109389/people/josie/NFDS_analysis/ballermix_holi11/outputs/formatted_freqs
FAI=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/NFDS/STAR.holi11.analysis.genome.fasta.fai

i=$(cat /gpfs/ts0/projects/Research_Project-T109423/people/bonnie/NFDS/STAR.holi11.analysis.genome.fasta.fai  | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1)

python /lustre/home/bf299/BallerMixPlus/BalLeRMix+_v1.py  --input $INPUT_DIR/${i}_APHP.derived.frq.derived --spect APHP_derived_spect.txt --o APHP_B2_${i}_non_win.out --usePhysPos --rec 2.12E-06 


