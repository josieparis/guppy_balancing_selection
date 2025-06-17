library(Rfast)
library(parallel)
library(data.table)


ghp <- fread('Documents/balancing_selection/vcftools_out/P.hwe', header = T)
head(ghp)
ghp_het <- strsplit(ghp$`OBS(HOM1/HET/HOM2)`, "/")
ghp_het_df <- as.data.frame(matrix(unlist(ghp_het), ncol=3, byrow=TRUE))
head(ghp_het_df)
colnames(ghp_het_df) <- c("HOM1", "HET", "HOM2")
ghp_het_df$HOM1 <- as.numeric(ghp_het_df$HOM1)
ghp_het_df$HET <- as.numeric(ghp_het_df$HET)
ghp_het_df$HOM2 <- as.numeric(ghp_het_df$HOM2)
ghp$total <- (ghp_het_df$HOM1 + ghp_het_df$HET + ghp_het_df$HOM2)
head(ghp)
ghp$HET_prop <- (ghp_het_df$HET/ghp$total)
hist(ghp$HET_prop)

chrs<-unique(ghp$CHR)
wind_size<-50000
#x="chr1"


#num of sites per window
winds<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-ghp[ghp$CHR == x,]
  winds1<-seq(0,max(tmp$POS),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$POS <= winds2[y] & tmp$POS >= winds1[y],]
    out<-data.frame(mean(tmp2$HET_prop))
    out$sum <- sum(tmp2$HET_prop != 0)
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'P'
    colnames(out)<-c('mean_het','num_het', 'chrom','window','BP1','BP2','comp')
    
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))

head(winds)
hist(winds$num_het/50000)

write.table(winds, "Documents/balancing_selection/vcftools_out/P_hwe_50k_winds.txt", sep = "\t", row.names=F, quote=F)


test <- subset(winds, chrom == "chr1")

chrom_plot <- ggplot(test, aes(x=(BP1+500), y=num_het))+
  geom_point(size = 0.5)+
#  geom_hline(yintercept = 100, linetype = "dashed", color = "red")+
#  scale_color_gradient(low = "pink", high = "darkorchid", limits = c(-3, 9))+
  #ylim(0,300)+
 # scale_fill_viridis_c(limits = c(-3, 9), oob = scales::squish, name = "Sum")+
  scale_x_continuous(name = "pos (MB)", labels=c(0,5,10,15,20,25,30,35,40,45), breaks=c(0,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000))+
  theme_bw()


chrom_plot

