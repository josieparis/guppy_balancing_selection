## Directory for ancestral allele liftover

This is how I lifted the AA alleles for the file `STAR_AA.final.fa`

To lift the AA, I used Emiliano Trucchi's ([emitruc](https://github.com/emitruc)) script (slightly modified).

NB a new (and recommended) version of this script is available at https://github.com/emitruc/genoloader/ 

- This script outputs the variants to be incorporated into the ALT column
- If they’re REF and AA they move to ALT, if they’re ALT and AA they also move to ALT.
- Code was added to count and output how many alleles (relative to the reference) been flipped
- The AA allele is now in the ALT column
- Also, need to record the AA allele when it is in fact REF.

#### Run the python script

NB needs python2 and scipy

`VCF=<dir>/holi11.SNP.maxmiss50_picta_wingei.vcf`

NB the vcf file has to be uncompressed otherwise you will get an error message  

`python2 ../scripts/writePolarizedVcf_toALT_replicateREF.py -f $VCF -p1 data/guppy.pop -p2 data/wingei.pop -p3 data/picta.pop`

(where p1 is the focal pop, p2 is the sister and p3 is the outgroup) - for each .pop file one individual per line

this will create a file called $VCF.pol.vcf with the ALT allele as AA, which you can then use for bcftools consensus to make a AA fasta:

`VCF_pol=vcf.pol.gz`

bgzip and tabix the pol.vcf:

```
bgzip ${VCF}_pol
tabix ${VCF}_pol.gz
```

set the reference genome:

`REF=STAR.chromosomes.release.fasta`

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

NB for the holi13 datatset (which had 4,532,226 SNPs, 4,298,080 had AA alleles, so 94%. This is what is in the STAR_AA.fasta)

