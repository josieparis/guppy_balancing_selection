### SNAPPER analyses

### VCF preparation:

`VCF=holi11_wingei_SNP.filtered.recode.vcf.gz`

1. Check for individuals with lowest amount of missing data:

`module load VCFtools`

`vcftools --vcf $VCF --missing-indv --out ${VCF}.missingness`

##### Choose three individuals with the lowest amount of missing data from each population. Save in file "SNAPP.inds"

2. Keep only inds and chrs we want

`SNAPP_popmap=SNAPP.inds`

`vcftools --vcf ${VCF} --keep ${SNAPP_popmap} --not-chr 12 --recode --recode-INFO-all --out SNAPP.inds`

3. Remove missing data and filter for mac

`vcftools --vcf SNAPP.inds.recode.vcf --max-missing 1.0 --min-mac 3 --recode --recode-INFO-all --out SNAPP.inds.nomiss.mac3`

4. Thin VCF (note that LD pruning here is inappropriate due to the inclusion of _P. wingei_)

`vcftools --vcf SNAPP.inds.nomiss.mac3.recode.vcf --thin 5000 --recode --recode-INFO-all --out SNAPP.inds.nomiss.mac3.thin5000`

##### NB different thinning thresholds were tested. 5k gave the highest nunber of useable SNPs for running beast 3 independent times (see below)

### Now ready to run SNAPPER

#### SNAPPER preparation

`snapp_VCF=SNAPP.inds.wing.bi.nomiss.mac3.thin5000.recode.vcf`

Using Michael Matschiner's ([mmatschiner](https://github.com/mmatschiner)) Ruby script: https://github.com/mmatschiner/snapp_prep 

```
$HOME/programs/matschiner/snapp_prep.rb \
-a SNAPPER -v ${snapp_VCF} \
-t poec_wing_samples.txt \
-c poec_wing_constraints.txt \
-s starting_tree_wing.nwk \
-l 1000000 -m 3000 -q 10000 \
-x poec_wing_snapper.xml \
-o poec_wing_snapper
```

#### Running snapper
```
$HOME/programs/beast/bin/beast \
-threads 8 \
poec_wing_snapper.xml
```
#### We ran 3 different analyses with different sets of SNPs
