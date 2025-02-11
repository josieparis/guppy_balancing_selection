

## SNAPP analysis
wd: `/lustre/home/jrp228/NERC/people/josie/NFDS_analysis/SNAPP`

### VCF
`VCF=/lustre/home/jrp228/startup/STAR_holi_snp_processing/holi11_vcf_files/holi11.SNP.maxmiss80_wingei.vcf`

### keep three individuals per pop (least missing data)
`vcftools --vcf ${VCF} --keep inds_to_keep.tsv --recode --recode-INFO-all --out holi11_3inds`

### remove missing data
`vcftools --vcf holi11_3inds.recode.vcf --max-missing 1 --recode --recode-INFO-all --out holi11_3inds_nomissing`

### calculate LD:
`vcf=holi11_3inds_nomissing`

`plink --vcf ${vcf}_IDs.vcf --out ${vcf}_pruned --indep-pairwise 50 5 0.2 --allow-extra-chr --double-id`

### convert to ped/map
`plink --bfile holi11_3inds_nomissing_NoLD --recode --out holi11_3inds_nomissing_NoLD --allow-extra-chr --double-id`

### Prune
`plink --ped holi11_3inds_nomissing.ped --map holi11_3inds_nomissing.map --extract holi11_3inds_nomissing_pruned.prune.in --out holi11_3inds_nomissing_pruned.vcf --recode vcf --allow-extra-chr --double-id`

### Set new VCF file:
`vcf=holi11_3inds_nomissing_pruned.vcf`

### Load Ruby
`module load Ruby/2.3.1`

### Run Micha's Ruby script:

https://github.com/mmatschiner/snapp_prep

`snapp_path=~/programs/snapp_prep`
