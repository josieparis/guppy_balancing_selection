library(data.table)
library(parallel)

B2_non_win<-read.table("Documents/balancing_selection/ballermix/non_win/B2/APHP_B2_ALL_non_win.out")
colnames(B2_non_win) <- c("physPos", "genPos", "CLR", "x_hat", "s_hat", "A_hat", "nSites","chr")
B2_non_win$log_s_hat <- log10(B2_non_win$s_hat)
summary(B2_non_win)

## summarise SNP values into windows
chrs<-unique(B2_non_win$chr)
wind_size<-50000

winds<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-B2_non_win[B2_non_win$chr == x,]
  winds1<-seq(0,max(tmp$physPos),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$physPos <= winds2[y] & tmp$physPos >= winds1[y],]
    if (nrow(tmp2) > 0) {
    max_row <- tmp2[which.max(tmp2$CLR), ]  # Find the row with the maximum CLR
    out <- data.frame(max_CLR = max_row$CLR, logS = max_row$log_s_hat)  # Create a data frame with max CLR and corresponding logS
    } else {
      out <- data.frame(max_CLR = NA, logS = NA)  # Handle empty windows
    }
    # Tidy
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'P'
    colnames(out)<-c('max_CLR','logS','chr','window','BP1','BP2','comp')
    
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))

write.table(winds, file = "Documents/balancing_selection/ballermix/non_win/B2/APHP_B2_ALL_non_win_windowised_100k.out", quote=F, row.names=F)

