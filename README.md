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
this means the AA allele is now in the ALT column
to run this script:

`module load Python/2.7.14-foss-2018a` (you need scipy) installed in your python environment

`VCF=

python ../scripts/writePolarizedVcf_toALT_and_count.py -f $VCF -p1 guppy.pop -p2 wingei.pop -p3 picta.pop

(where p1 is the focal pop, p2 is the sister and p3 is the outgroup)

this will create a file called $VCF.pol.vcf with the ALT allele as AA, which you can then use for bcftools consensus to make a AA fasta:

`bcftools consensus -f $REF picta_wingei_combined_newAA.vcf.gz -o STAR_AA.fa`

`samtools faidx STAR_AA.fa`

and then lift:

`zcat holi_13_4532226_final.vcf.gz | ~/programs/vcftools/src/perl/fill-aa -a STAR_AA.fasta | bgzip -c > holi13_pictawingei_AA.vcf.gz`




