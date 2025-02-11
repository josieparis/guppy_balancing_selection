### Script to run through PicMin picking of ballermix results
# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures


## load libs
lib <- c("data.table","ggplot2","qvalue","Rfast","pbmcapply","poolr","dplyr", "tidyr", "tibble")

sapply(lib,library,character.only=T)
n_cores = 6


#' @title orderStatsPValues
#' @param p_list the vector of p-values from your genome scans
#' @importFrom "stats" "pbeta"
orderStatsPValues <- function(p_list){
  ## This function returns a list of the p-values for each of marginal p-values
  # sort the list of $p$-values
  p_sort <- sort(p_list)[2:length(p_list)]
  # calculate the number of species minus 1
  n = length(p_sort)
  # get a vector of the 'a' parameters for each of the marginal distributions
  the_as = 2:(length(p_sort) + 1)
  # get a vector of the 'b' parameters for each of the marginal distributions
  the_bs = n+2-the_as
  # calculate the p-values for each of the marginals and return
  return(pbeta(p_sort, the_as, the_bs))
}

#' @title generatePvalCorrMat
#' @param N the number of comparisons being made, including the first pval that is discarded by picmin
#' @param geneN the number of genes being simulated
#' @param seed set a seed integer for reproducible random correlation matrices
#' @importFrom "Rfast" "rowSort"
#' @importFrom "stats" "pbeta"
generatePvalCorrMat <- function(N,geneN=10000,seedN = 1000){
  set.seed(seedN)
  # Draw a matrix where each row is a 'gene' and each column is a replicate species/lineage
  ee <- matrix(nrow = geneN,ncol = N)
  for (i in 1:N){
    temp <- runif(geneN)
    ee[,i] <- EmpiricalPs(temp)
  }
  # Sort these row-wise from smallest to largest
  ee_sorted <- Rfast::rowSort(ee)
  # Column-wise, estimate the pval from the beta distribution
  beta_ps = matrix(nrow = geneN,ncol = N - 1)
  for(i in 1:ncol(beta_ps)){
    j = i + 1
    beta_ps[,i] <- pbeta(ee_sorted[,j],j,(N - j + 1))
  }
  # Calculate correlation matrix among cols
  cor_ord_unif <- cor(beta_ps)
  cor_ord_unif
}

generatePvalCorrMat_BUG <- function(N,geneN=10000,seed = 1000){
  set.seed(seed)
  ee <- array (NA,c(geneN,N))
  for (i in 1:n){
    temp <- runif(10000)
    ee[,i] <- EmpiricalPs(temp)
  }
  emp_p_null_dat_unif <- ee
  
  for (kk in 1:nrow(emp_p_null_dat_unif)){
    ind1 <- order(emp_p_null_dat_unif[kk,])
    emp_p_null_dat_unif [kk,] <- emp_p_null_dat_unif[kk,ind1]
  }
  cor_ord_unif <- cor (emp_p_null_dat_unif[,2:N])
  cor_ord_unif
}

#' @title PicMin
#' @param N the number of comparisons being made, including the first pval that is discarded by picmin
#' @param geneN the number of genes being simulated
#' @param seed set a seed integer for reproducible random correlation matrices
#' @param numReps the number of replicate draws to perform when building the empirical distributing for calculating the Tippett p-value 
#' @importFrom "poolr" "empirical"
generateTippettNull <- function(N,geneN = 10000,seedN = 1000,numReps = 100000){
  tmp_cor_mat = generatePvalCorrMat(N = N,geneN = geneN,seed = seedN)
  poolr:::empirical(R = tmp_cor_mat,
                    method = 'tippett',
                    side = 1,
                    size = numReps)
}


#' @title PicMin
#' @param pList the vector of p-values from your genome scans
#' @param correlationMatrix the correlation matrix under the null hypothesis
#' @param numReps the number of replicate draws to perform when building the empirical distributing for calculating the Tippett p-value 
#' @importFrom "poolr" "tippett"
PicMin <- function(pList, 
                   correlationMatrix, 
                   numReps = 100000,
                   fixConfig = NULL){
  # Calculate the p-value for the order statistics
  ord_stats_p_values <- orderStatsPValues(pList)
  # Apply the Tippett Correction
  if(is.null(fixConfig)){
    p_value <- tippett(ord_stats_p_values,
                       adjust = "empirical",
                       R = correlationMatrix,
                       side = 1, size = numReps)$p
    return(list(p=p_value,
                config_est=which.min(ord_stats_p_values)+1))
  } else if(length(fixConfig) > 1) {
    # Modified version of PicMin to handle looking at a subset of configs
    p_value <- tippett(ord_stats_p_values[fixConfig - 1], 
                       adjust = "empirical", 
                       R = correlationMatrix[fixConfig - 1,fixConfig - 1], 
                       side = 1, size = numReps)$p
    return(list(p = p_value,
                config_est = fixConfig[which.min(ord_stats_p_values[fixConfig - 1])]))
  } else {
    # Modified version of PicMin to handle looking at a specific config
    p_value <- ord_stats_p_values[fixConfig - 1]
    return(list(p=p_value,
                config_est=fixConfig))
  }
}

#' @title PicMin_bulk
#' @param pvals the vector of p-values from your genome scans across all groups being tested
#' @param groups the vector of group values, for e.g. orthogroup, gene, or window_id
#' @param nulls a list of empirical null distributions calculated with generateTippettNull(). Should be a separate empirical null for each set size being tested. The list must be named according to the set size. If only one set size is being tested, i.e. all groups have N replicates, a single vector can be provided.
#' @import "data.table"
#' @importFrom "qvalue" "empPvals"
PicMin_bulk = function(pvals,
                       groups,
                       nulls){
  
  # Make a new data.table
  test_dd = data.table(test_group = groups,pval = pvals)
  
  # First assign the a's and b's, also remove the lowest pval from each
  ordered = test_dd[,.(a_param = rank(.SD$pval,ties.method = 'random'),
                       # b_param = (nrow(.SD) - rank(.SD$pval) + 1),
                       testN = nrow(.SD),
                       pval),by = test_group][a_param != 1,]
  ordered$b_param = ordered$testN - ordered$a_param + 1
  # ordered$b_param = ordered$testN - ordered$a_param
  
  # Calculate the beta param distribution and retain the minimum and do correction
  ordered$beta_p = pbeta(ordered$pval,ordered$a_param,ordered$b_param)
  ordered_DS = unique(ordered[,.(minP_DS = 1 - (1 - min(.SD$beta_p))^nrow(.SD),
                                 testN,
                                 config_est = .SD$a_param[which.min(.SD$beta_p)]),by = test_group])
  
  # Tippett correction of DS-corrected pval using null
  testN_vals = sort(unique(ordered_DS$testN))
  if(length(testN_vals) > 1){
    rbindlist(lapply(testN_vals,function(i){
      tmp = ordered_DS[testN == i,]
      tmp$picmin_p = qvalue::empPvals(-tmp$minP_DS,-nulls[[ as.character(i) ]])
      tmp[,.(test_group,picmin_p,config_est,testN)]
    }))
  } else {
    ordered_DS$picmin_p = qvalue::empPvals(-ordered_DS$minP_DS,-nulls)
    ordered_DS[,.(test_group,picmin_p,config_est,testN)]
  }
  
}

#' @title Generate data under the null hypothesis
#' @param adaptation_screen The threshold used to determine adaptation
#' @param a the 'a' parameter of a beta distribution of p-values for the false null
#' @param b the 'b' parameter of a beta distribution of p-values for the false null
#' @param n the number of species in the test
#' @param genes the number of genes in the genome use to calculate empirical p-values
#' @importFrom "stats" "rbeta"
GenerateNullData <- function(adaptation_screen, a, b, n, genes){
  temp <- c( rbeta(1,a,b),replicate(n-1, sample(genes,1)/genes) ) 
  while (sum(temp<adaptation_screen)==0){
    temp <- c( rbeta(1,a,b),replicate(n-1, sample(genes,1)/genes) ) 
  }
  return(temp)
}

#' @title Calculate empirical p-values from a vector of summary statistics
#' @param vector_of_values A set of summary statistics that you want to convert to empirical p-values
#' @param large_i_small_p Do you want large values to have small p-values (e.g. Fst)?
EmpiricalPs <- function( vector_of_values, large_i_small_p = FALSE ){
  if  (large_i_small_p==TRUE){
    rank(vector_of_values * -1,na.last = "keep")/sum(is.na (vector_of_values) == F)
  }
  else{
    rank(vector_of_values,na.last = "keep")/sum(is.na (vector_of_values) == F)
  }
}

#' @title Generate N random genes to use for making a null distribution of picmin results under randomness
#' @param Ngenes The number of genes that are being tested in the observed data
#' @param Nsims The number of simulated sets of genes
#' @param Ngroups The number of species/populations/groups per gene
#' @import "data.table"
GenerateRandomGenes <- function(Ngenes,
                                Nsims = 1e5,
                                Ngroups){
  
  # Make a matrix and populate with random draws of empirical pvalues of the order we expect
  random_matrix = matrix(ncol = Nsims,
                         nrow = Ngroups)
  for(i in 1:nrow(random_matrix)){
    random_matrix[i,] = sample((1:Ngenes)/Ngenes,Nsims,replace = T)
  }
  
  # Melt this and return for PicMin
  out = reshape2::melt(random_matrix)
  colnames(out) = c("group","gene","pval")
  data.table(out)
}

##here is where you start putting in your data###

glist <-c("chr1","chr2", "chr3", "chr4", "chr5", "chr6", "chr7","chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23")


for(i in glist){
  chr <- i
  input_name <- paste("Documents/balancing_selection/ballermix/non_win/B2/picmin/ALL_pops_B2_non_win_",chr,".csv", sep = "")
  

  dd <- read.csv(file = input_name)
 
   colnames(dd) <- c("physPos", "genPos","CLR","x_hat","s_hat","A_hat", "nSites","chr","pop")
  dd$ID <- paste(dd$chr, dd$physPos, sep = "_")


# Count all the IDs
ID_count <- table(dd$ID)
to_run <- names(ID_count[ID_count == max(ID_count)])
dd <- dd[dd$ID %in% to_run,]

# Emp P matrix
#pop to comp.x

emp_mat <- matrix(ncol=length(unique(dd$comp.x)),nrow=length(to_run))
for(i in 1:length(unique(dd$comp.x))){
  
  pop_tmp <- dd[dd$comp.x == unique(dd$comp.x)[i],]
  pop_tmp <- pop_tmp[order(pop_tmp$ID),]
  emp_mat[,i] <- empPvals(pop_tmp$max_CLR,pop_tmp$max_CLR)
  
}
colnames(emp_mat) <- unique(dd$comp.x)
rownames(emp_mat) <- pop_tmp$ID
nrow(emp_mat)

emp_mat <- as.data.frame(emp_mat)

reshaped_emp_mat <- emp_mat %>%
  rownames_to_column(var = "RowNames") %>%
  pivot_longer(cols = -RowNames, names_to = "Variable", values_to = "Value")

reshaped_emp_mat <- as.data.frame(reshaped_emp_mat)

head(reshaped_emp_mat)

########starting new workflow from here######
# Assuming here you have 11 populations with balancing selection results
testN = 11
windowN = nrow(emp_mat)
repN = 10000000 # How many reps do you want for tippett correction, assume 1 million

# Set up the Tippett Nulls
tippett_nulls = generateTippettNull(N = testN,
                                    geneN = windowN,
                                    numReps = repN,
                                    seedN = 123)

# Set up data for input
# This should be some data.table or something where one column is the pvals in a window
# and the other column says which window

picmin_input <- data.frame(pval = reshaped_emp_mat$Value, window_id = reshaped_emp_mat$RowNames)


# Run PicMin
# I made this wrapper for vectorising picmin
# Only inputs needed are a vector of pvals, a vector saying what group each pval is in
# and a tippett null distribution
picmin_res = PicMin_bulk(pvals = picmin_input$pval,
                         groups = picmin_input$window_id,
                         nulls = tippett_nulls)



#summary(picmin_res)
# FDR correct and tidy
#picmin_res$picmin_fdr = p.adjust(picmin_res$picmin_p_adj,method = 'fdr')

picmin_res$picmin_fdr = p.adjust(picmin_res$picmin_p,method = 'fdr')

picmin_res$pos <- sapply(strsplit(picmin_res$test_group,"_"),'[[',2)


results_name <- paste("Documents/balancing_selection/ballermix/non_win/B2/picmin/Picmin_res",chr,"txt", sep = ".")
write.table(picmin_res, file = results_name, sep = "\t", row.names=F, quote=F)

}
