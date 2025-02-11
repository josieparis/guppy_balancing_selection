### Liftover for Chr20 and Chr12
# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures


setwd("/Users/josie/Dropbox/Sussex_Guppies/Analyses/NFDS_analysis/vcf_files")


source("/Users/josie/Dropbox/Sussex_Guppies/Analyses/NFDS_analysis/R_functions/chr20_scaf94_liftover.R")
source("/Users/josie/Dropbox/Sussex_Guppies/Analyses/NFDS_analysis/R_functions/chr12_liftover.R")

lib<-c("data.table","vcfR","tidyverse","parallel","ggplot2")
lapply(lib,library,character.only=T)

###############
#### LG12 ####
###############

vcf_file <- "holi11_chr12_shapeit_beagle.vcf.gz"

vcf<-read.vcfR(vcf_file)
meta<-data.frame(vcf@fix)

# Run the liftover function for the new chr12 STAR positions:
meta[meta$CHROM == "chr12","POS"]<-update_STAR(as.integer(meta[meta$CHROM == "chr12","POS"]))
# Replace metadata
vcf@fix<-as.matrix(meta)

# Write new chr12 VCF
write.vcf(file = "holi11_chr12_phased.update.vcf.gz", vcf)

## need to sort the vcf file afterwards using:
# bcftools sort LG12.AA.update.vcf.gz -Oz -o LG12.AA.update.sorted.vcf.gz
# Then rename to poolseq_variant_filtered_update_STAR.vcf.gz

#########################
#### LG20 / scaff 94 ####
########################

## Example usage
# merge_scaf94_chr20(dd[dd$chr == "chr20","BP"],scaf="chr20")
# merge_scaf94_chr20(dd[dd$chr == "000094F_0","BP"],scaf="94")

vcf_file <- "holi11_chr20_scaf94.vcf.gz"

vcf<-read.vcfR(vcf_file)
meta<-data.frame(vcf@fix)

# Transform the BP locations
fix<-data.frame(vcf@fix)
fix$POS<-as.character(fix$POS)
chr20<-as.integer(fix[fix$CHROM == "chr20","POS"])
scaf94<-as.integer(fix[fix$CHROM == "000094F_0","POS"])

# Turn into new
chr20_new<-merge_scaf94_chr20(chr20,"chr20")
scaf94_new<-merge_scaf94_chr20(scaf94,"94")


# Replace
fix[fix$CHROM == "chr20","POS"]<-as.character(as.integer(chr20_new))
fix[fix$CHROM == "000094F_0","POS"]<-as.character(as.integer(scaf94_new))
fix[fix$CHROM == "000094F_0","CHROM"]<-"chr20"

# Sort
vcf@fix<-as.matrix(fix[order(c(chr20_new,scaf94_new)),])
vcf@gt<-vcf@gt[order(c(chr20_new,scaf94_new)),]

# Write the new vcf
write.vcf(vcf,
          "holi11.chr20_scaf94_phased.update.vcf.gz")

system("rm -f tmp.vcf")





