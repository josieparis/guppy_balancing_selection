### Rscript for popGenome fst, pi and tajD ###
#install.packages("ff",lib='~/bin/R/library',repos='https://www.stats.bris.ac.uk/R/')
#install.packages("data.table",lib='~/bin/R/library',repos='https://www.stats.bris.ac.uk/R/')
#install.packages("dplyr",lib='~/bin/R/library',repos='https://www.stats.bris.ac.uk/R/')
#install.packages("PopGenome",lib='~/bin/R/library',repos='https://www.stats.bris.ac.uk/R/')

library("ff", lib.loc='~/bin/R/library')
library("data.table", lib.loc='~/bin/R/library')
library("PopGenome", lib.loc='~/bin/R/library')

##UPDATE WITH FAI FOR WHATEVER VCF YOURE USING
#Import fai for genome data
chr_length<-read.table("/gpfs/ts0/projects/Research_Project-T109423/STAR_holi_snp_processing/vcf_files_for_analysis/STAR.holi11_analysis.reference.fasta.fai",header=F)

####### Loop over chromosomes ###########
chr_function<-function(x){

#chr name as defined from fasta.fai
tid=as.character(chr_length[x,]$V1)


#chr length
topos = chr_length[x,]$V2


#Read in VCF
vcf_file <- readVCF('/gpfs/ts0/projects/Research_Project-T109423/STAR_holi_snp_processing/vcf_files_for_analysis/holi11.rearranged.vcf.gz', tid=tid, frompos = 1, topos = topos, numcols = 1000000, include.unknown=TRUE)
### assign pops ###
ghp<-c('GH1','GH10','GH11','GH12','GH13','GH14','GH15','GH17','GH18','GH19','GH2','GH20','GH3','GH4','GH5','GH6','GH7','GH8','GH9')
glp<-c('GL1','GL10','GL11','GL12','GL13','GL14','GL15','GL16','GL17','GL18','GL19','GL2','GL20','GL5','GL6','GL7','GL8','GL9')
tuhp<-c('TUHP_F1','TUHP_F10','TUHP_F2','TUHP_F3','TUHP_F4','TUHP_F5','TUHP_F6','TUHP_F7','TUHP_F8','TUHP_F9','TUHP_M1','TUHP_M10','TUHP_M2','TUHP_M3','TUHP_M4','TUHP_M5','TUHP_M6','TUHP_M7','TUHP_M8')
tulp<-c('TULP_F1','TULP_F10','TULP_F2','TULP_F3','TULP_F4','TULP_F5','TULP_F6','TULP_F7','TULP_F8','TULP_F9','TULP_M1','TULP_M10','TULP_M2','TULP_M3','TULP_M4','TULP_M5','TULP_M6','TULP_M7','TULP_M8','TULP_M9')
aphp <-c('APHP_F1','APHP_F10', 'APHP_F2', 'APHP_F3', 'APHP_F4', 'APHP_F5', 'APHP_F6',' APHP_F7', 'APHP_F8', 'APHP_F9', 'APHP_M1', 'APHP_M10', 'APHP_M2', 'APHP_M3', 'APHP_M5', 'APHP_M6', 'APHP_M7', 'APHP_M8', 'APHP_M9')
aplp <-c('APLP_F1', 'APLP_F10', 'APLP_F3', 'APLP_F4', 'APLP_F5', 'APLP_F6', 'APLP_F7', 'APLP_F8', 'APLP_F9', 'APLP_M1', 'APLP_M10', 'APLP_M2','APLP_M3', 'APLP_M4', 'APLP_M5', 'APLP_M6', 'APLP_M7', 'APLP_M8')
echp <-c('ECHP_F1', 'ECHP_F10', 'ECHP_F2', 'ECHP_F3', 'ECHP_F4', 'ECHP_F5', 'ECHP_F6', 'ECHP_F7', 'ECHP_F8', 'ECHP_F9', 'ECHP_M1', 'ECHP_M10', 'ECHP_M2','ECHP_M3','ECHP_M4','ECHP_M5','ECHP_M6','ECHP_M7','ECHP_M9')
eclp <-c('ECLP_F1','ECLP_F10','ECLP_F2','ECLP_F3','ECLP_F4','ECLP_F5','ECLP_F7','ECLP_F8', 'ECLP_F9','ECLP_M1','ECLP_M10', 'ECLP_M2','ECLP_M3','ECLP_M4','ECLP_M5','ECLP_M6','ECLP_M7','ECLP_M8', 'ECLP_M9')
mhp <- c('LM1','LM10','LM11','LM12','LM13','LM14','LM15','LM16','LM17','LM18','LM19','LM2','LM20','LM3','LM4','LM5','LM6','LM7','LM8','LM9')
mlp <-c('UM1','UM10','UM11','UM12','UM13','UM14','UM15','UM16','UM17','UM18','UM2','UM3','UM5','UM6','UM7','UM8','UM9')
p <- c('Paria_F12','Paria_F16','Paria_F2','Paria_F2_RNA','Paria_F8','Paria_M12','Paria_M14','Paria_M2','Paria_M5','Paria_M9')


vcf_file <- set.populations(vcf_file, list(ghp,glp,tuhp,tulp,aphp,aplp,echp,eclp,mhp,mlp,p), diploid=TRUE)

### create windows ###
slide_vcf_file_75kb <- sliding.window.transform(vcf_file, 50000, 50000, type=2)

### do pop stats
slide_vcf_file_75kb <- F_ST.stats(slide_vcf_file_75kb, mode="nucleotide")
slide_vcf_file_75kb <- neutrality.stats(slide_vcf_file_75kb, FAST=TRUE)
slide_vcf_file_75kb <- diversity.stats.between(slide_vcf_file_75kb, nucleotide.mode=TRUE)
slide_vcf_file_75kb <- detail.stats(slide_vcf_file_75kb,site.spectrum=TRUE)


 calc_FST<-function(fst_list){
         tmp<-data.frame(GHP_GLP= fst_list[,1],
 		GHP_TUHP= fst_list[,2],
 		GHP_TULP= fst_list[,3],
 		GHP_APHP= fst_list[,4],
 		GHP_APLP= fst_list[,5],
 		GHP_ECHP= fst_list[,6],
 		GHP_ECLP= fst_list[,7],
  		GHP_MHP=fst_list[,8], 
  		GHP_MLP=fst_list[,9],
  		GHP_P=fst_list[,10],
 		GLP_TUHP= fst_list[,11],
 		GLP_TULP= fst_list[,12],
 		GLP_APHP= fst_list[,13],
 		GLP_APLP= fst_list[,14],
 		GLP_ECHP= fst_list[,15],
 		GLP_ECLP= fst_list[,16],
  		GLP_MHP=fst_list[,17],
 		GLP_MLP=fst_list[,18],
 		GLP_P=fst_list[,19],
  		TUHP_TULP= fst_list[,20],
 		TUHP_APHP= fst_list[,21],
		TUHP_APLP= fst_list[,22],
		TUHP_ECHP= fst_list[,23],
		TUHP_ECLP= fst_list[,24],
		TUHP_MHP=fst_list[,25],
		TUHP_MLP=fst_list[,26],
		TUHP_P=fst_list[,27],
		TULP_APHP= fst_list[,28],
		TULP_APLP= fst_list[,29],
		TULP_ECHP= fst_list[,30],
		TULP_ECLP= fst_list[,31],
		TULP_MHP=fst_list[,32],
		TULP_MLP=fst_list[,33],
		TULP_P=fst_list[,34],
  		APHP_APLP=fst_list[,35],
  		APHP_ECHP=fst_list[,36],
  		APHP_ECLP=fst_list[,37],
  		APHP_MHP=fst_list[,38],                
    		APHP_MLP=fst_list[,39],
		APHP_P=fst_list[,40],
  		APLP_ECHP=fst_list[,41],
  		APLP_ECLP=fst_list[,42],
  		APLP_MHP=fst_list[,43],
    		APLP_MLP=fst_list[,44],
    		APLP_P=fst_list[,45],
		ECHP_ECLP= fst_list[,46],
		ECHP_MHP=fst_list[,47],
    		ECHP_MLP=fst_list[,48],
     		ECHP_P=fst_list[,49],
    		ECLP_MHP=fst_list[,50],
     		ECLP_MLP=fst_list[,51],
     		ECLP_P=fst_list[,52],
     		MHP_MLP=fst_list[,53],
       		MHP_P=fst_list[,54],
       		MLP_P=fst_list[,55])

      if(length(tmp) == 0){
         return('NA')
         } else {
         tmp['chrom']<-rep(tid,nrow(tmp))
         tmp['window']<-1:nrow(tmp)
         tmp['window_start']<- (tmp$window-1)*50000
         tmp['window_end']<- tmp$window*50000
         tmp<-tmp[,c(56:59,1:55)]
         }
         }
 FST_mat<-data.frame(calc_FST(t(slide_vcf_file_75kb@nuc.F_ST.pairwise)))
outfile_fst <-paste0("/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/NFDS/popgenome/",tid,".50kb.fst.popgenome.out")
 write.table(FST_mat, outfile_fst, sep = "\t", row.names=F, quote = F)

 calc_PI<-function(pi_list){
         tmp<-data.frame(GHP = pi_list[,1]/50000,
                 GLP = pi_list[,2]/50000,
                 TUHP = pi_list[,3]/50000,
                 TULP = pi_list[,4]/50000,
		 APHP = pi_list[,5]/50000,
		 APLP = pi_list[,6]/50000,
		 ECHP = pi_list[,7]/50000,
		 ECLP = pi_list[,8]/50000,
		 MHP = pi_list[,9]/50000,
		MLP = pi_list[,10]/50000,
		P = pi_list[,11]/50000)
         if(length(tmp) == 0){
         return('NA')
         } else {
         tmp['chrom']<-rep(tid,nrow(tmp))
         tmp['window']<-1:nrow(tmp)
         tmp['window_start']<- (tmp$window-1)*50000
         tmp['window_end']<- tmp$window*50000
         tmp<-tmp[,c(12:15,1:11)]
         }
         }
 PI_mat<-data.frame(calc_PI(slide_vcf_file_75kb@nuc.diversity.within))
 outfile_pi <-paste0("/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/NFDS/popgenome/",tid,".50kb.pi.popgenome.out")
 write.table(PI_mat, outfile_pi, sep = "\t", row.names=F, quote = F)

 calc_TD<-function(td_list){

         tmp<-data.frame(GHP= td_list[,1],
                 GLP = td_list[,2],
                 TUHP = td_list[,3],
		 TULP = td_list[,4],
		 APHP = td_list[,5],
                 APLP = td_list[,6],
                 ECHP = td_list[,7],
                 ECLP = td_list[,8],
                 MHP = td_list[,9],
                MLP = td_list[,10],
                P = td_list[,11])

         if(length(tmp) == 0){
         return('NA')
         } else {
         tmp['chrom']<-rep(tid,nrow(tmp))
         tmp['window']<-1:nrow(tmp)
         tmp['window_start']<- (tmp$window-1)*50000
         tmp['window_end']<- tmp$window*50000
         tmp<-tmp[,c(12:15,1:11)]
         }
         }
 TD_mat<-data.frame(calc_TD(slide_vcf_file_75kb@Tajima.D))
 outfile_td <-paste0("/gpfs/ts0/projects/Research_Project-T109423/people/bonnie/NFDS/popgenome/",tid,".50kb.td.popgenome.out")
 write.table(TD_mat, outfile_td, sep = "\t", row.names=F, quote = F)

}

x_vector<-seq(1,nrow(chr_length))
lapply(x_vector,chr_function)
