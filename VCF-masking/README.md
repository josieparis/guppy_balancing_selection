### VCF filtering prior to scans of balancing selection

#### NB Matches Supplementary Methods S3

`VCF=holi11.SNP.maxmiss50.maf0.03_AA.tags`

I want to mask for the following:
1) Low mappability
2) High repeat content
3) Excess low or high coverage variants

### Mappability:

The positive mappability mask is STAR.mappability_mask.positive.fasta.gz

`map_bed=STAR_mappability_mask.positive.bed`

Mask with bedtools intersect:

```
bedtools intersect -a ${VCF}.vcf -b STAR_mappability_mask.positive.bed -wa -f 1.0 > ${VCF}.mapmask.vcf
```

### Repeats:

`repeat_bed=STAR_repeats_10kb_percentage.bed`

Plot distribution:

<img width="848" alt="Screenshot 2023-06-07 at 19 50 21" src="https://github.com/josieparis/NFDS/assets/38511308/1d327849-960a-4397-b8d9-b5a2149c635e">

Choose 15% repeats used as threshold for high repeat content

Filter:

`cat STAR_repeats_10kb_percentage.bed | awk '{if ($5>=15) {print}}' > STAR_repeats_over15%.bed`

`STAR_repeats_over15%.bed | wc -l` # 4580 windows

The bedtools reverse intersect on this bed file:

```
bedtools intersect -a ${VCF}.mapmask.vcf -b STAR_repeats_over15%.bed -v > ${VCF}.mapmask.repeatmask.vcf
```

### Coverage:

Assessed from the VCF file, removing SNPs with coverage lower than and higher than a threshold:

```
vcftools --vcf ${VCF}.mapmask.repeatmask.vcf --site-mean-depth --out ${VCF}.mapmask.repeatmask
```

Plot distribution:

<img width="844" alt="Screenshot 2023-06-12 at 17 48 35" src="https://github.com/josieparis/NFDS/assets/38511308/a149633e-f312-47ab-a40e-1e3767fcaa4a">


7 used as a min low depth coverage and 20 used as a mean high coverage filter:

```
vcftools --vcf holi11.SNP.maxmiss50.maf0.03_AA.tags_map_repeat_filtered.vcf --minDP 7 --max-meanDP 20 --minGQ 30 --max-missing 0.5 --recode --out ${VCF}.mapmask.repeatmask.coveragemask
```


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
