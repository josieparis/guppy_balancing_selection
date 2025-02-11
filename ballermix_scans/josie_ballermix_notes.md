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
