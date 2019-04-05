###This is a script to make small files out of the big 'combined.sure.files' that now each contain all data for each chromosome
##these will be easier to handle in the downstream analysis


setwd("~/mydata/git/SuRE_shared/Joris/SNP_project/all_snp_analysis_wilcox/")


##get the file names for all 46 files (23 chromosomes but for both strands)
file.vector <- paste('~/mydata/git_data/Downsampled_SuRE_counts_180802/combined.sure.snp.files.180802/',list.files(path = '~/mydata/git_data/Downsampled_SuRE_counts_180802/combined.sure.snp.files.180802/'),sep='')
file.vector <- file.vector[2:47]

file.vector

###so i want to know the snp positions and select based on that

#install.packages("foreach")
library(foreach)

#install.packages("doMC")
library(doMC)
registerDoMC(1)  #change the 2 to your number of CPU cores  



foreach(k=1:23) %dopar% {
  print(paste('starting file number ', k))
  
##first combine the files from the two strands of the same chromosome  
sure.df.1 <- readRDS(file.vector[k])
sure.df.2 <- readRDS(file.vector[k+23])


sure.df <- rbind(sure.df.1,sure.df.2)

##order by SNP position
sure.df.ordered <- sure.df[order(sure.df$SNPabspos),]


##when you split the file by a certain number of SNPs, you want the file to contain the data from all elements containing those SNP
##when the files are ordered on absolute position of the SNPs a new SNP starts every time when you find a duplicated position. 


new.snps <- which(!duplicated(sure.df.ordered$SNPabspos))
tail(new.snps)

##make blocks of 1000 snps.

new.snp.positions <- new.snps[seq(from=1,to=length(new.snps),by=1000)]
snp.last.positions <- c(new.snp.positions-1,nrow(sure.df.ordered))
blocks <-snp.last.positions 

###now save those blocks of 1000 SNPs

for(i in 1:(length(blocks)-1)){
  print(i)
  saveRDS(sure.df.ordered[(blocks[i]+1):blocks[i+1],],paste('~/mydata/git_data/Downsampled_SuRE_counts_180802/combined.sure.snp.files.180802/combined.sure.snp.files.small/',sub('.2018.*','',sub('.*snp.combined.minus','sure.snp.combined',file.vector[k])),'.',Sys.Date(),'.',i,'.downsampled.rds',sep=''))
  
}

}
