## GitHub repo containing data and scripts pertaining to the article:

### Repeated signatures of balancing selection in small and large populations of guppies (_Poecilia reticulata_)

#### Authors
Josephine R Paris, James R Whiting, Joan Ferrer Obiol, Kimberly A Hughes, Bonnie A Fraser

#### Available on bioRxiv: xxx

### Archived under the Zenodo DOI: xxx

---------

### VCF files used for analysis:

MASTER directory for SNP calling:
`/lustre/home/jrp228/startup/STAR_holi_snp_processing`

MASTER directory with holi11 vcf files:
`/lustre/home/jrp228/startup/STAR_holi_snp_processing/holi11_vcf_files`

This includes the data on APHP, APLP, MHP, MLP, Paria, GHP, GLP, TUHP, TULP, ECHP, ECLP

VCF files for different analyses:

- PCA: holi11.SNP.maxmiss80.vcf_IDs.vcf (649,984 variants)
- fineStructure: holi11.SNP.maxmiss80.vcf.gz (phased version - see below)
- relate: holi11.SNP.maxmiss50.vcf.gz (no maf filter - phased version - see below)
- popgen stats (including Taj D): holi11.SNP.maxmiss50.maf0.03.recode.vcf.gz (3% maf filter) - 3,502,712
- Ballermix: holi11.SNP.maxmiss50_picta_wingei.maf0.03.vcf (3% maf filter, has integrated AA alleles from picta / wingei data (see below)) - 3,378,040
- Baypass: holi11.SNP.maxmiss50.maf0.03.recode.vcf.gz - 3,502,712


Phased vcfs:
`/lustre/home/jrp228/startup/STAR_holi_snp_processing/holi11_vcf_files/phased_vcfs`

Phased vcf for fineStructure (no maf filter): `holi11_vcf_files/phased_vcfs/holi11_maxmiss80_nomaf`

Phased vcf for relate (no maf filter): `holi11_vcf_files/phased_vcfs/holi11_maxmiss50_nomaf`

--------- 

Notes (to be deleted or ordered in another file):

## Betascan:
NB: Info on how to run Betascan is detailed in main lab book - will copy these lines over for clarity  
I have decided to rerun it on the holi11 vcf file. I also want to mask the variants, as per the ballermix paper (masking below)

To use the B2 statistic we need information on the ancestral allele (so we can calculate the derived allele frequency).
This comes from the picta and wingei data. The picta and wingei data are with holi11 in this vcf file:

`holi11.SNP.maxmiss50_picta_wingei.vcf.gz` 

and for the balancing selection scans, I'm using the maf 3% filtered one (see above)

NB: this file has **biallelic** variants (so cannot do multiallelic bal sel scans)

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

NB obviously scaffolds where the stats aren't calculated will give a warning message.

This also appends all the header lines from each of the data frames, we want to keep one (for the header line) and remove all the others:

`grep "^chrom" holi11_allchrs.50kb.td.out | head -n 1 > header`

`grep -v "^chrom" holi11_allchrs.50kb.td.out > tmp`

`echo $'chrom\twindow\twindow_start\twindow_end\tAPHP\tAPLP\tECHP\tECLP\tGHP\tGLP\tMHP\tMLP\tTUHP\tTULP\tPARIA' | cat - tmp > holi11_allchrs.50kb.td.out`

`rm tmp batch_inputs.txt header`

---------------











