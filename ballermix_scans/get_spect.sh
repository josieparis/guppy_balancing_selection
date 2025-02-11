#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 
#SBATCH -A Research_Project-T109423
#SBATCH -p pq
#SBATCH --job-name=ballermix
#SBATCH --error=ballermix.err.txt
#SBATCH --output=ballermix.out.txt
#SBATCH --export=All
#SBATCH -D .

module load SciPy-bundle/2020.03-foss-2020a-Python-3.8.2
INPUT_DIR=/gpfs/ts0/projects/Research_Project-T109389/people/josie/NFDS_analysis/ballermix_holi11/outputs/formatted_freqs

#loop for all pops
pop_array=(APHP APLP ECHP ECLP GHP GLP MHP MLP P TUHP TULP)

for POP in "${pop_array[@]}"
do


#concat all pop files without header line
awk 'FNR>1' $INPUT_DIR/chr*_${POP}.derived.frq.derived > ${POP}_concat_derived_frq.txt
awk 'BEGIN {print "phsPos\tgenPos\tx\tn"} {print}' ${POP}_concat_derived_frq.txt > temp.txt && mv temp.txt ${POP}_concat_derived_frq.txt

#get spect
python /lustre/home/bf299/BallerMixPlus/BalLeRMix+_v1.py --getSpect --input ${POP}_concat_derived_frq.txt --spect ${POP}_derived_spect.txt

done
