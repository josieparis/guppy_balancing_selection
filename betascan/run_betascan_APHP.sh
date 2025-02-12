#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16 # nodes=number of nodes required. ppn=number of processors per node
#SBATCH -A Research_Project-T109423
#SBATCH -p pq
#SBATCH --job-name=betascan
#SBATCH --error=betascan.err.txt
#SBATCH --output=betascan.out.txt
#SBATCH --export=All
#SBATCH -D .
#SBATCH --array=1-25%25

module load SciPy-bundle/2020.03-foss-2020a-Python-3.8.2
INPUT_DIR=/gpfs/ts0/projects/Research_Project-T109389/people/josie/NFDS_analysis/ballermix_holi11/outputs/formatted_freqs
FAI=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/NFDS/STAR.holi11.analysis.genome.fasta.fai


#then make an array per chr
i=$(cat /gpfs/ts0/projects/Research_Project-T109423/people/bonnie/NFDS/STAR.holi11.analysis.genome.fasta.fai  | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1)

awk 'NR > 1 {print $1, $3, $4}' $INPUT_DIR/${i}_APHP.derived.frq.derived  > APHP_out/${i}_APHP_betascan_in.txt

python /lustre/home/bf299/BetaScan-master/BetaScan.py -i APHP_out/${i}_APHP_betascan_in.txt -o APHP_out/${i}_APHP_betascan_out
