# NFDS-paper
Thoughts, scripts, analysis  for NFDS balancing selection paper

I am making a new README file to keep all the NFDS information together (out of main lab book). Can also share this with Bonnie.

---------

### VCF files for analysis:

MASTER directory for SNP calling:
`/lustre/home/jrp228/startup/STAR_holi_snp_processing`

MASTER directory with holi11 vcf files:
`/lustre/home/jrp228/startup/STAR_holi_snp_processing/holi11_vcf_files`

This includes the data on APHP, APLP, MHP, MLP, Paria, GHP, GLP, TUHP, TULP, ECHP, ECLP

Vcf files for different analyses:

- PCA: holi11.SNP.maxmiss80.vcf_IDs.vcf (649,984 variants)
- fineStructure: holi11.SNP.maxmiss80.vcf.gz (phased version - see below)
- relate: holi11.SNP.maxmiss50.vcf.gz (no maf filter - phased version - see below)
- popgen stats (including Taj D): holi11.SNP.maxmiss50.maf0.03.recode.vcf.gz (3% maf filter)
- Ballermix: holi11.SNP.maxmiss50_picta_wingei.maf0.03.vcf (3% maf filter)
- Baypass: holi11.SNP.maxmiss50.maf0.03.recode.vcf.gz

Phased vcfs:
`/lustre/home/jrp228/startup/STAR_holi_snp_processing/holi11_vcf_files/phased_vcfs`

Phased vcf for fineStructure (no maf filter): `holi11_vcf_files/phased_vcfs/holi11_maxmiss80_nomaf`

Phased vcf for relate (no maf filter): `holi11_vcf_files/phased_vcfs/holi11_maxmiss50_nomaf`

-----------

## Betascan:
NB: Info on how to run Betascan is detailed in main lab book - will copy these lines over for clarity  
I have decided to rerun it on the holi11 vcf file. I also want to mask the variants, as per the ballermix paper (masking below)

To use the B2 statistic we need information on the ancestral allele (so we can calculate the derived allele frequency).
This comes from the picta and wingei data. The picta and wingei data are with holi11 in this vcf file:

`holi11.SNP.maxmiss50_picta_wingei.vcf.gz` 

and for the balancing selection scans, I'm using the maf 3% filtered one (see above)

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

`python2 ../scripts/writePolarizedVcf_toALT_replicateREF.py -f $VCF -p1 data/guppy.pop -p2 data/wingei.pop -p3 data/picta.pop`

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

`VCF=/lustre/home/jrp228/startup/STAR_holi_snp_processing/holi11_vcf_files/holi11.SNP.maxmiss50_picta_wingei.vcf`

`python scripts/writePolarizedVcf_toALT_replicateREF.py -f $VCF -p1 data/guppy.pop -p2 data/wingei.pop -p3 data/picta.pop` 

Number of flipped ancestral alleles:  1,358,374 (so out of 16,213,978, the command flipped 8.4% of alleles)

load BCFtools:

`module load BCFtools/1.9-foss-2018a`

bgzip the pol.vcf:

`bgzip ${VCF}_pol`

tabix:

`tabix ${VCF}_pol.gz`

set the STAR reference genome:

`REF=~/startup/STAR/STAR.chromosomes.release.fasta`

Run bcftools consensus

`bcftools consensus -f $REF ${VCF}_pol.gz -o ${VCF}_AA.fa`

index:
`samtools faidx $VCF}_AA.fa`

Load the Perl libraries required for fill-aa:

`export PERL5LIB=~/programs/vcftools/src/perl/`

load samtools:

module load SAMtools/1.9-foss-2018b

and then lift:

NB here the $VCF is the original VCF you need to add the AA fields into (not the polarised one)

`zcat $VCF | ~/programs/vcftools/src/perl/fill-aa -a ${VCF}_AA.fasta | bgzip -c > ${VCF}_AA.vcf.gz`

Now fill the tags with bcftools:

`bcftools +fill-tags $VCF_AA.vcf.gz | bgzip -c > ${VCF}_AA.tags.vcf.gz`

NB this will fill the tags for the main VCF file, but if you then split by pop (as below in step 2), the AA tag will disappear again, so you need to add it back in below

#### Step 2:
Now you need to generate an AA VCF file per population:

`pops=(APHP APLP MHP MLP P GHP GLP TUHP TULP ECHP ECLP)`

this is done in a submission script, e.g.:

`vcftools --gzvcf holi11.SNP.maxmiss50.maf0.03_AA.tags.vcf.gz --keep ./popmaps/APHP.popmap --recode --recode-INFO-all --out APHP.AA`

then bgzip and tabix all the outputs:

`parallel bgzip {} ::: *.vcf`

`parallel tabix {} ::: *.vcf.gz`

fill the tags:

`for pop in ${pops[@]}; do bcftools +fill-tags ${pop}.AA.recode.vcf.gz | bgzip -c > ${pop}.AA.tags.vcf.gz; done`

Now run the vcftools freq option on each of these (asking for the derived allele) using the script 01_calc_derived_freqs.sh:

```
## Load modules
module purge
module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0

MASTER=/lustre/home/jrp228/NERC/people/josie/NFDS_analysis/ballermix_holi11
vcf_files=/lustre/home/jrp228/startup/STAR_holi_snp_processing/holi11_vcf_files/population_vcfs
freqs=$MASTER/outputs/01_raw_derived_freqs

## Set up the chroms
#chrs=$(awk '{print $1}' ${reference}.fai | sed "${SLURM_ARRAY_TASK_ID}q;d")

## Set up the pops
#pops=(APHP APLP MHP MLP P GHP GLP TUHP TULP ECHP ECLP)
pops=(MLP P GHP GLP TUHP TULP ECHP ECLP)

for pop in ${pops[@]}
do
vcftools --gzvcf $vcf_files/${pop}.AA.tags.vcf.gz --chr chr${SLURM_ARRAY_TASK_ID} --freq --derived --out ${freqs}/chr${SLURM_ARRAY_TASK_ID}_${pop}.derived
done
```

Concurrently, run easySFS.py to assess the best sample size to use for the spect SFS (column N_CHR in the freqs output).

Note that you have to run these per chromosome as ballermix takes a per chromosome input. Also each chromosome has it's own averaged recombination rate (derived from Jim's Heredity paper)

Note from this post: https://github.com/bioXiaoheng/BalLeRMix/issues/5 that you can use one spectrum file for the whole genome (you don't need one per chromosome), but you need to use the input derived alleles for each chromosome


---------------
## Tajima's D PicMin analysis

- run PicMin on the whole genome, not per chromosome
- Our popgenome script obviously outputs for each chromosome separately, so first you need to concatenate:

`reference=/lustre/home/jrp228/startup/STAR/STAR.chromosomes.release.fasta`

`chr_array=($(awk '{print $1}' ${reference}.fai))`

`for i in "${chr_array[@]}"; do echo "$i.50kb.td.popgenome.out" >> batch_inputs.txt; done`

`xargs -i cat '{}' < batch_inputs.txt > holi11_allchrs.50kb.td.out`

NB obviously scaffolds where the stats aren't calculated will give a warning message




---------------
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







