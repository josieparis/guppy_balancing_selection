## Directory for ancestral allele liftover

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

`cd  ../STAR_holi_snp_processing/holi11_vcf_files/ancestral_allele_vcfs`

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
`samtools faidx ${VCF}_AA.fa`

Load the Perl libraries required for fill-aa:

`export PERL5LIB=~/programs/vcftools/src/perl/`

load samtools:

module load SAMtools/1.9-foss-2018b

and then lift:

NB here the $VCF is the original VCF you need to add the AA fields into (not the polarised one)

`zcat $VCF | ~/programs/vcftools/src/perl/fill-aa -a ${VCF}_AA.fasta | bgzip -c > ${VCF}_AA.vcf.gz`

Now fill the tags with bcftools:

`bcftools +fill-tags $VCF_AA.vcf.gz | bgzip -c > ${VCF}_AA.tags.vcf.gz`

NB this will fill the tags for the main VCF file, but if you then split by pop (as below in step xx), the AA tag will disappear again, so you need to add it back in below


#### Step 4:
Now you need to generate an AA VCF file per population:

`pops=(APHP APLP MHP MLP P GHP GLP TUHP TULP ECHP ECLP)`

this is done in a submission script, e.g.:

`vcftools --gzvcf holi11.SNP.filtered.AA.tags.vcf --keep ./popmaps/APHP.popmap --recode --recode-INFO-all --out APHP.AA`

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


---------------------------------------- 
## Scaffolds analysis


Also will extract information for all scaffolds > 500,000 bp:

```
000032F_0.2     707072
000077F_0       2478161
000083F_0.3     766951
000094F_0       1797025
000095F_0       1773698
000104F_0       1540729
000111F_0       1427647
000113F_0       1367121
000117F_0       1249272
000119F_0       1209297
000122F_0       1109139
000126F_0       1036673
000135F_0       847025
000140F_0       746010
000149F_0       630695
000150F_0       611721
000151F_0       627290
000152F_0       603738
000153F_0       586083
000154F_0       571063
000155F_0       564286

```

Analysis is inside the folder: `/lustre/home/jrp228/startup/STAR_holi_snp_processing/holi11_vcf_files/scaffold_extraction`

Number of flipped ancestral alleles:  16,239 (out of 222,532)

AA sites filled  .. 79098
AA bases filled  .. 79098


