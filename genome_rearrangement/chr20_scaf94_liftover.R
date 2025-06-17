################################################################################################
# Based on HiC information, merge scaffold 94 with chr20 

merge_scaf94_chr20<-function(x,scaf=NULL){
  
  ## Example usage
  # merge_scaf94_chr20(dd[dd$chr == "chr20","BP"],scaf="chr20")
  # merge_scaf94_chr20(dd[dd$chr == "000094F_0","BP"],scaf="94")
  
  merged_out<-sapply(x,function(bp){
    
    # Is it on SCF_204?
    if(scaf=="chr20"){
      if(bp > 836423 & bp <= 3164071){
        bp_out<- 3164071 - bp + 836423 + 1797025
      } else if (bp < 836423) {
        bp_out<-bp
      } else {
        bp_out<-bp+1797025
      }
    } else if (scaf=="94"){
      bp_out<-1797025 - bp + 836423
    }
    return(bp_out)
  })
  return(as.integer(merged_out))
}

