update_STAR <- function(bp_vec=NULL){
  vec_out <- sapply(bp_vec,function(x){  
    # Invert I and II...
    if(x <= 3258871){
      bp_out <- 3258871 - x
    } 
    # Move Contig 6...
    else if (x >= 5331818 & x <= 7997830){
      bp_out <- x + 433864 + 10000
    } else if (x >= 8007831 & x <= 8441694){
      bp_out <- x - 2666013 - 10000
    }
    # Invert Contig XI
    else if (x >= 20209082 & x <= 24252560){
      bp_out <- 24252560 - x + 20209082
    }
    # Swap Contig XIII and XIV
    else if (x >= 24511363 & x <= 25631364){
      bp_out <- x + 963876 + 10000
    } else if (x >= 25641365 & x <= 26605240){
      bp_out <- x - 1120002 - 10000
    } else {
      bp_out <- x
    }
    return(bp_out)
  })
  return(vec_out)
}


