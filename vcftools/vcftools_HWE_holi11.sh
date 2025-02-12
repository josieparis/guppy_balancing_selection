#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -A Research_Project-T109423
#SBATCH -p pq
#SBATCH --job-name=vcf_HWE
#SBATCH --error=vcf_HWE.err.txt
#SBATCH --output=vcf_HWE.out.txt
#SBATCH --export=All
#SBATCH -D .

module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0
in_vcf=/gpfs/ts0/projects/Research_Project-T109423/STAR_holi_snp_processing/vcf_files_for_analysis/holi11.rearranged.vcf.gz
out_dir=/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/NFDS/vcftools
pop_files=/gpfs/ts0/projects/Research_Project-T109423/STAR_holi_snp_processing/holi11_vcf_files/popmaps

### getting allele frequencies
vcftools --gzvcf $in_vcf \
--keep $pop_files/APHP.popmap \
--hardy \
--out $out_dir/APHP

vcftools --gzvcf $in_vcf \
--keep $pop_files/APLP.popmap \
--hardy \
--out $out_dir/APLP

vcftools --gzvcf $in_vcf \
--keep $pop_files/ECHP.popmap \
--hardy \
--out $out_dir/ECHP

vcftools --gzvcf $in_vcf \
--keep $pop_files/ECLP.popmap \
--hardy \
--out $out_dir/ECLP

vcftools --gzvcf $in_vcf \
--keep $pop_files/GHP.popmap \
--hardy \
--out $out_dir/GHP

vcftools --gzvcf $in_vcf \
--keep $pop_files/GLP.popmap \
--hardy \
--out $out_dir/GLP

vcftools --gzvcf $in_vcf \
--keep $pop_files/MHP.popmap \
--hardy \
--out $out_dir/MHP

vcftools --gzvcf $in_vcf \
--keep $pop_files/MLP.popmap \
--hardy \
--out $out_dir/MLP

vcftools --gzvcf $in_vcf \
--keep $pop_files/P.popmap \
--hardy \
--out $out_dir/P

vcftools --gzvcf $in_vcf \
--keep $pop_files/TUHP.popmap \
--hardy \
--out $out_dir/TUHP

vcftools --gzvcf $in_vcf \
--keep $pop_files/TULP.popmap \
--hardy \
--out $out_dir/TULP

