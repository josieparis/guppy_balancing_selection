# NFDS
Thoughts, scripts etc for NFDS bal sel paper

I am making a new README file to keep all the NFDS information together (out of main lab book). Can also share this with Bonnie.

### Vcf files for analysis:

MASTER directory for SNP calling:
`/lustre/home/jrp228/startup/STAR_holi_snp_processing`

MASTER directory with holi11 vcf files:
`/lustre/home/jrp228/startup/STAR_holi_snp_processing/holi11_vcf_files`

This includes the data on APHP, APLP, MHP, MLP, Paria, GHP, GLP, TUHP, TULP, ECHP, ECLP

Vcf files for different analyses:

- PCA: holi11.SNP.maxmiss80.vcf_IDs.vcf (649,984 variants)
- fineStructure: holi11.SNP.maxmiss80.vcf.gz (phased version - see below)
- relate: holi11.SNP.maxmiss50.vcf.gz (no maf filter - phased version - see below)
- popgen stats: holi12_19563701.nomaf.vcf.gz (no maf filter)

Phased vcfs:
`/lustre/home/jrp228/startup/STAR_holi_snp_processing/holi12_vcf_files/phased_vcfs`

Phased vcf for fineStructure (no maf filter): `holi11_vcf_files/phased_vcfs/holi11_maxmiss80_nomaf` 
Phased vcf for relate (no maf filter): `holi11_vcf_files/phased_vcfs/holi11_maxmiss50_nomaf`

## Betascan:
NB: Info on how to run Betascan is detailed in main lab book - will copy these lines over for clarity  
I have decided to rerun it on the holi11 vcf file. I also want to mask the variants, as per the ballermix paper (masking below)

---
## Step 1 
To use the B2 statistic we need information on the ancestral allele
This comes from the picta and wingei data. The picta and wingei data are with holi11 in this vcf file:

`holi11.SNP.maxmiss50_picta_wingei.vcf.gz` (16,213,978 variants)
NB: this file has **biallelic** variants (so cannot do multiallelic bal sel scans)

This is how I lifted the AA alleles for the file `/lustre/home/jrp228/startup/STAR/STAR_AA.final.fa` (this was done with the holi13_picta_wingei vcf file so maximum number of samples (across all our studies) was included - this was intended to be a "generic" AA fasta for our other studies.

This was how it was done:

To lift the AA, I modified emiliano's script so that it outputs the variants to be incorporated into the ALT column. so if they’re REF and AA they move to ALT, if they’re ALT and AA they also move to ALT. Then I added some code to count and output how many have been flipped
this means the AA allele is now in the ALT column. Also, need to record the AA allele when it is in fact REF.

#### Step 1:

`module load Python/2.7.15-foss-2018b` (you also need scipy) installed in your python environment. I did this with `pip install scipt --user`

`VCF=/lustre/home/jrp228/startup/STAR_holi_snp_processing/holi11_vcf_files/holi11.SNP.maxmiss50_picta_wingei.vcf.gz`

NB the vcf file has to be unzipped otherwise you will get the error message:

```
Traceback (most recent call last):
  File "scripts/writePolarizedVcf_toALT_and_count.py", line 156, in <module>
    index_sample = header.index(sample)
NameError: name 'header' is not defined
```

`cd  ... /STAR_holi_snp_processing/holi11_vcf_files/ancestral_allele_vcfs`

python2 ../scripts/writePolarizedVcf_toALT_replicateREF.py -f $VCF -p1 data/guppy.pop -p2 data/wingei.pop -p3 data/picta.pop

(where p1 is the focal pop, p2 is the sister and p3 is the outgroup) - for each .pop file one individual per line

this will create a file called $VCF.pol.vcf with the ALT allele as AA, which you can then use for bcftools consensus to make a AA fasta:

`VCF_pol=vcf.pol.gz`

bgzip the pol.vcf:

`bgzip ${VCF}_pol`

tabix:

`tabix ${VCF}_pol.gz`

set the STAR reference genome:

`REF=~/startup/STAR/STAR.chromosomes.release.fasta`

load BCFtools:

`module load BCFtools/1.9-foss-2018a`

#### Run bcftools consensus

`bcftools consensus -f $REF ${VCF}_pol.gz -o ${VCF}_AA.fa`

index:
`samtools faidx $VCF}_AA.fa`

Load the Perl libraries required for fill-aa:

`export PERL5LIB=~/programs/vcftools/src/perl/`

and then lift:

NB here the $VCF is the original VCF you need to add the AA fields too (not the polarised one)

`zcat $VCF | ~/programs/vcftools/src/perl/fill-aa -a ${VCF}_AA.fasta | bgzip -c > ${VCF}_AA.vcf.gz`

Now fill the tags with bcftools:

`bcftools +fill-tags $VCF_AA.vcf.gz | bgzip -c > ${VCF}_AA.tags.vcf.gz`

NB for the holi13 datatset (which had 4,532,226 SNPs, 4,298,080 had AA alleles, so 94%) (this is what is in the STAR_AA.fasta)

### Applying this to holi11:

```
VCF=/lustre/home/jrp228/startup/STAR_holi_snp_processing/holi11_vcf_files/holi11.SNP.maxmiss50_picta_wingei.vcf

python scripts/writePolarizedVcf_toALT_and_count.py -f $VCF -p1 data/guppy.pop -p2 data/wingei.pop -p3 data/picta.pop 

Number of flipped ancestral alleles:  1,358,374

```





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

`







