## GitHub repo containing data and scripts pertaining to the article:

### Repeated signatures of balancing selection in small and large populations of guppies (_Poecilia reticulata_)

#### Authors
Josephine R Paris, James R Whiting, Joan Ferrer Obiol, Kimberly A Hughes, Bonnie A Fraser

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

NOTES on chromosome / scaffold rearrangements:
- For LG12, you need to reorder the contigs. May also need to add in scaffold 149 after Contig XIII:
<img width="1072" alt="Screenshot 2023-06-07 at 15 36 29" src="https://github.com/josieparis/NFDS/assets/38511308/7577a976-45a8-4114-81e5-38af51a50591">

- For LG20, we may need to add scaffold 94 at the beginning
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


#### Step 3: 
I now want to mask the VCF file: holi11.SNP.maxmiss50.maf0.03_AA.tags.vcf.gz

I will mask for the following:
1. mappability mask
2. high / low coverage regions
3. repeat regions

Before I started, I included the scaffolds (see top of document) which had also been ancestral allele flipped. They are here:

```/lustre/home/jrp228/startup/STAR_holi_snp_processing/holi11_vcf_files/scaffold_extraction/scaffolds_holi13.maxmiss50.maf0.03.AA_tags.vcf.gz```

I concatenated the scaffolds with the file above, so the file is now called "concatenated.vcf.gz"

3,457,138 SNPs

##### Mappability:

the mappability mask is here: `/lustre/home/jrp228/startup/mappability`

there's a positive fasta (i.e. nucleotides as P pass mappablity, those with N do not): `STAR.mappability_mask.positive.fasta`

the bed file that made this mask is here also: `kmer_intersect_filtered.bed`

Mask with bedtools intersect:

```
bedtools intersect -a concatenated.vcf.gz -b kmer_intersect_filtered.bed -wa -f 1.0 > tmp.vcf
```

3,447,879 SNPs

This VCF then needs to be intersected with high repeat regions removed
##### Repeats:

Are here: 
`/lustre/home/jrp228/startup/STAR/STAR_repeats_10kb_percentage.bed`

I copied this file to the holi11_vcf folder

<img width="848" alt="Screenshot 2023-06-07 at 19 50 21" src="https://github.com/josieparis/NFDS/assets/38511308/1d327849-960a-4397-b8d9-b5a2149c635e">

15% repeats used as threshold for high repeat content ...

I filtered the all_repeats_count_10kb.bed file to those windows with more than 15% repeats 

`cat STAR_repeats_10kb_percentage.bed | awk '{if ($5>=15) {print}}' > STAR_repeats_over15%.bed`

STAR_repeats_over15%.bed | wc -l # 4580 windows

The bedtools reverse intersect on this bed file:
```bedtools intersect -a tmp2.vcf -b STAR_repeat_mask.bed -v > tmp3.vcf```

3,246,395

NB:
get the header:

zgrep "^#" concatenated.vcf.gz > vcf_header

paste it onto the tmp3 file and rename the VCF

cat vcf_header tmp3.vcf > tmp4.vcf


##### Coverage:
Can be assessed from the VCF file, removing SNPs with coverage lower than and higher than a threshold:

```vcftools --vcf tmp4.vcf --site-mean-depth --out tmp4```

visualise in geom_density:

<img width="844" alt="Screenshot 2023-06-12 at 17 48 35" src="https://github.com/josieparis/NFDS/assets/38511308/a149633e-f312-47ab-a40e-1e3767fcaa4a">


5 used as a min low depth coverage and 20 used as a mean high coverage filter:

```
vcftools --vcf holi11.SNP.maxmiss50.maf0.03_AA.tags_map_repeat_filtered.vcf --minDP 7 --max-meanDP 16 --minGQ 30 --max-missing 0.5 --recode --out holi11.SNP.maxmiss50.maf0.03_AA.tags_map_repeat_cov_filtered
```

Outputting VCF file...
After filtering, kept 3234526 out of a possible 3246395 Sites


### Summary of filtering
```
Starting VCF   3,502,712  holi11.SNP.maxmiss50.maf0.03.recode.vcf.gz
With AA field 3,378,040 holi11.SNP.maxmiss50.maf0.03_AA.vcf.gz (-124,672)
With the scaffolds 3,457,138 concatenated.vcf.gz
Mappability mask  3,447,879 SNPs tmp.vcf
Repeat mask 3,246,395 tmp3.vcf
Coverage mask 3,234,526 tmp4.vcf
```

clean up:
```
mv holi11.SNP.maxmiss50.maf0.03_AA.tags_map_repeat_cov_filtered.recode.vcf holi11.SNP.filtered.AA.vcf
bgzip holi11.SNP.filtered.AA.vcf
tabix holi11.SNP.filtered.AA.vcf.gz
rm tmp.vcf
rm tmp2.vcf
rm vcf_header
rm holi11.SNP.maxmiss50.maf0.03_AA.tags_map_repeat_filtered.vcf
```

Finally, redo the AA liftover and filling of AA tags:
(in the ancestral_allele_vcf folder):
```
cat holi11.SNP.filtered.AA.vcf | ~/programs/vcftools/src/perl/fill-aa -a STAR_AA.final.fa | bgzip -c > holi11.SNP.filtered_actuallywithAA.vcf.gz
bcftools +fill-tags holi11.SNP.filtered_actuallywithAA.vcf.gz | bgzip -c > holi11.SNP.filtered_actuallywithAA_tags.vcf.gz
```

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

Concurrently, run easySFS.py to assess the best sample size to use for the spect SFS (column N_CHR in the freqs output).

Note that you have to run these per chromosome as ballermix takes a per chromosome input. Also each chromosome has it's own averaged recombination rate (derived from Jim's Heredity paper)

Note from this post: https://github.com/bioXiaoheng/BalLeRMix/issues/5 that you can use one spectrum file for the whole genome (you don't need one per chromosome), but you need to use the input derived alleles for each chromosome



---------------------------------------- 
## Scaffolds analysis

Analysis is inside the folder: `/lustre/home/jrp228/startup/STAR_holi_snp_processing/holi11_vcf_files/scaffold_extraction`

Number of flipped ancestral alleles:  16,239 (out of 222,532)

AA sites filled  .. 79098
AA bases filled  .. 79098




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

## Ballermix scans August 2023

### Main vcf file (contigs rearranged and with ancestral alleles):

`/lustre/home/jrp228/startup/STAR_holi_snp_processing/vcf_files_for_analysis/holi11.rearranged.vcf.gz`

### ballermix has to be run on a per chromosome / population basis. The population-specific VCF files are here:

`/lustre/home/jrp228/startup/STAR_holi_snp_processing/vcf_files_for_analysis/pop_specific`

NB: I checked that these VCF files all have the scaffolds and the AA column for calculating the derived allele frequency 

### Then we need to calculate the derived allele frequency across each chromosome in each population:

(There is a script for this, called 01_calc_derived_freqs.sh which is here:

'/lustre/home/jrp228/NERC/people/josie/NFDS_analysis/ballermix_holi11/scripts'

Run this for all chromosomes and scaffolds 

### Then we need to format the derived allele counts so that they work for Ballermix. You can do this in the freqs folder, giving the output to another folder:

`for i in chr*frq; do awk -F "\t|:" '(NR>1) && ($8!='0') && ($3=='2') {OFS="\t"; print$2,$8*$4,$4}' OFMT="%.0f" $i > ../formatted_freqs/$i.derived; done`

This code removes the monomorphic sites (column 8 is 0), makes sure the sites are biallelic (column 2 is 2) and then prints the position column, the allele frequency * the sample size and then the sample size. 

Add the NA column for the physical position:
`for i in *derived; do awk '{OFS="\t"; {print $1,"NA",$2,$3}}' $i > $i.derived2; done`

`rename .derived.derived2 .derived *derived2`

Add the header line:
`for i in *derived; do echo -e 'phsPos\tgenPos\tx\tn' >> $i ; done`

make the last line the first line 
sed -i '1h;1d;$!H;$!d;G' *.derived

### Start using Ballermix

First, we need to generate the spectrum file for each chromosome, for each pop:

`for i in *derived; do python ~/programs/BalLeRMix_2.2/software/BalLeRMix_v2.5.py -i $i --getSpect --spect ../spect_files/spect.$i.DAF.txt; done`

rename:

`rename .out.frq.derived.DAF.txt .DAF.txt *txt`

#### Now the issue with the sample size still exists ... 
For now try with highest sample size:
APHP - 26
APLP - 24
ECHP - 34
ECLP - 36
GHP - 36
GLP - 36
MHP - 38
MLP - 34
P - 18 
TUHP - 32
TULP - 24

Generate these files as, e.g.: (for each pop)

`for i in chr*_TULP*derived; do awk '{ if ($4 == 24) { print } }' $i > ../highest_sample_size_freqs/$i.samplesize; done`










