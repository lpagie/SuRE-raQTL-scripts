##this script was used to perform the wilcoxon test for each SNP


setwd("~/mydata/git/SuRE_shared/Joris/SNP_project/all_snp_analysis_wilcox/")

#install.packages("foreach")
library(foreach)

#install.packages("doMC")
library(doMC)
registerDoMC(20)  

##refer to the small fiels generated of the combine SuRE-count tables
file.vector <- paste('~/mydata/git_data/Downsampled_SuRE_counts_180802/combined.sure.snp.files.180802/combined.sure.snp.files.small/',list.files(path = '~/mydata/git_data/Downsampled_SuRE_counts_180802/combined.sure.snp.files.180802/combined.sure.snp.files.small/'),sep='')
head(file.vector)


##loop through all files

foreach(k=1:length(file.vector)) %dopar% {
print(paste('starting file number ', k))
 
sure.df <- readRDS(file.vector[k])

head(sure.df)
sum(is.na(sure.df$sum.norm.K562.cDNA))

##elements will be grouped by their unique SNP position
snp.pos.vector <- unique(sure.df$SNPabspos)  


results.df <- data.frame(matrix(ncol=21, nrow=length(snp.pos.vector)))
colnames(results.df) <- c('chr','SNP_ID','SNPabspos','ref.element.count','alt.element.count','ref.element.plus.strand.count','alt.element.plus.strand.count','ref.element.minus.strand.count','alt.element.minus.strand.count','k562.ref.mean','k562.ref.median','k562.alt.mean','k562.alt.median','k562.wilcox.p.value','k562.wilcox.p.value.random','hepg2.ref.mean','hepg2.ref.median','hepg2.alt.mean','hepg2.alt.median','hepg2.wilcox.p.value','hepg2.wilcox.p.value.random')


##loop through all SNPs
for(i in 1: length(snp.pos.vector)){
  

  snp.idx <- which(sure.df $SNPabspos==snp.pos.vector[i])
  

  ref <- which(sure.df [snp.idx,'SNPvarInf']=="0") ##these are the elements containing the ref allele 
  alt <- which(sure.df [snp.idx,'SNPvarInf']=="1") ##these are the elements containing hte alt allele
  
  ##generate a random assignment of ref and alt of the same size as the real identities
  ref.random <- sample(c(ref,alt),length(ref))
  alt.random <- c(ref,alt)[!c(ref,alt)%in%ref.random]
  
  ##construct a results dataframe
  results.df[i,'SNP_ID'] <- sure.df [snp.idx[1],'SNP_ID']
  results.df[i,'SNPabspos'] <- sure.df [snp.idx[1],'SNPabspos']
  results.df[i,'ref.element.count'] <- length(ref) ##number of elements containing the ref allele
  results.df[i,'alt.element.count'] <- length(alt) ##number of elements containing the alt allele
  
  ##also save strand specific counts of elements for later checks
  results.df[i,'ref.element.plus.strand.count'] <- length(which(sure.df [snp.idx,'SNPvarInf']=="0" & sure.df [snp.idx,'strand']=="+"))
  results.df[i,'alt.element.plus.strand.count'] <- length(which(sure.df [snp.idx,'SNPvarInf']=="1" & sure.df [snp.idx,'strand']=="+"))
  results.df[i,'ref.element.minus.strand.count'] <- length(which(sure.df [snp.idx,'SNPvarInf']=="0" & sure.df [snp.idx,'strand']=="-"))
  results.df[i,'alt.element.minus.strand.count'] <- length(which(sure.df [snp.idx,'SNPvarInf']=="1" & sure.df [snp.idx,'strand']=="-"))
  
  
  
  results.df[i,'k562.ref.mean'] <- mean(sure.df [snp.idx[ref],'ipcr.norm.sum.K562.cDNA'])
  results.df[i,'k562.alt.mean'] <- mean(sure.df [snp.idx[alt],'ipcr.norm.sum.K562.cDNA'])
  results.df[i,'k562.ref.median'] <- median(sure.df [snp.idx[ref],'ipcr.norm.sum.K562.cDNA'])
  results.df[i,'k562.alt.median'] <- median(sure.df [snp.idx[alt],'ipcr.norm.sum.K562.cDNA'])
  
  ##perform wilcox but only if you have data for at least 1 ref element and 1 alt element 
  if(sum(is.finite(sure.df [snp.idx[ref],'ipcr.norm.sum.K562.cDNA']))>1 & sum(is.finite(sure.df [snp.idx[alt],'ipcr.norm.sum.K562.cDNA']))>1){
  results.df[i,'k562.wilcox.p.value'] <- wilcox.test(sure.df [snp.idx[ref],'ipcr.norm.sum.K562.cDNA'],sure.df [snp.idx[alt],'ipcr.norm.sum.K562.cDNA'])$p.value
  } else{}
  
  ##perform wilxoc again but using the randomly assigned alleles
  if(sum(is.finite(sure.df [snp.idx[ref.random],'ipcr.norm.sum.K562.cDNA']))>1 & sum(is.finite(sure.df [snp.idx[alt.random],'ipcr.norm.sum.K562.cDNA']))>1){
  results.df[i,'k562.wilcox.p.value.random'] <- wilcox.test(sure.df [snp.idx[ref.random],'ipcr.norm.sum.K562.cDNA'],sure.df [snp.idx[alt.random],'ipcr.norm.sum.K562.cDNA'])$p.value
  } else{}
    
  ##repeat what you did for the K562 data, now for the HepG2 data
  
  results.df[i,'hepg2.ref.mean'] <- mean(sure.df [snp.idx[ref],'ipcr.norm.sum.HEPG2.cDNA'])
  results.df[i,'hepg2.alt.mean'] <- mean(sure.df [snp.idx[alt],'ipcr.norm.sum.HEPG2.cDNA'])
  results.df[i,'hepg2.ref.median'] <- median(sure.df [snp.idx[ref],'ipcr.norm.sum.HEPG2.cDNA'])
  results.df[i,'hepg2.alt.median'] <- median(sure.df [snp.idx[alt],'ipcr.norm.sum.HEPG2.cDNA'])
    

  
  if(sum(is.finite(sure.df [snp.idx[ref],'ipcr.norm.sum.HEPG2.cDNA']))>1 & sum(is.finite(sure.df [snp.idx[alt],'ipcr.norm.sum.HEPG2.cDNA']))>1){
  results.df[i,'hepg2.wilcox.p.value'] <- wilcox.test(sure.df [snp.idx[ref],'ipcr.norm.sum.HEPG2.cDNA'],sure.df [snp.idx[alt],'ipcr.norm.sum.HEPG2.cDNA'])$p.value
   } else{}
  
  if(sum(is.finite(sure.df [snp.idx[ref.random],'ipcr.norm.sum.HEPG2.cDNA']))>1 & sum(is.finite(sure.df [snp.idx[alt.random],'ipcr.norm.sum.HEPG2.cDNA']))>1){
  results.df[i,'hepg2.wilcox.p.value.random'] <- wilcox.test(sure.df [snp.idx[ref.random],'ipcr.norm.sum.HEPG2.cDNA'],sure.df [snp.idx[alt.random],'ipcr.norm.sum.HEPG2.cDNA'])$p.value
  } else{}
  
  
}


results.df$chr <- paste('chr',sub('.2018.*','',sub('.*chr', '',file.vector[k])), sep='')
head(results.df)
saveRDS(results.df,paste('~/mydata/git_data/Downsampled_SuRE_counts_180802/wilcox.180806/',sub('.*combined.', '',file.vector[k]), sep=''))

}



###combine the files, add some columns, subset on coverage and save in entirety

file.names <- paste('~/mydata/git_data/Downsampled_SuRE_counts_180802/wilcox.180806/',list.files(path = "~/mydata/git_data/Downsampled_SuRE_counts_180802/wilcox.180806/"),sep='')

sure.snp.df <- readRDS(file.names[1])
for(i in 1:length(file.names)){
print(i)
temp.df <- readRDS(file.names[i])
sure.snp.df <- rbind(sure.snp.df,temp.df)
}




head(sure.snp.df)



##add columns for min/max coverage of elements and min/max expression
sure.snp.df$min.element <- apply(sure.snp.df[,c('ref.element.count','alt.element.count')],1,min)
sure.snp.df$max.element <- apply(sure.snp.df[,c('ref.element.count','alt.element.count')],1,max)
sure.snp.df$max.k562.expression <- apply(sure.snp.df[,c('k562.ref.mean','k562.alt.mean')],1,max)
sure.snp.df$max.hepg2.expression <- apply(sure.snp.df[,c('hepg2.ref.mean','hepg2.alt.mean')],1,max)

###save entire dataset, although i will later work with a subset

saveRDS(sure.snp.df, file='sure.snp.df.strand.nonspecific.exact.wilcox.downsampled.JvA.180806.rds')


###i will set a minimal and maximum threshold for coverage.
##these cutoffs were determined emperically. at lowecoverage the data becomes noisy, at very high coverage we get lots of contribution from amplified regions or repeats which somehow are all mapped to one part of the reference genome
covered.sure.snp.df <- sure.snp.df[which(sure.snp.df$min.element>4 & sure.snp.df$max.element<1000),]

###

covered.sure.snp.df <- covered.sure.snp.df[order(covered.sure.snp.df$k562.wilcox.p.value),]

##save the reduced file; this is the file used throughout the manuscript
saveRDS(covered.sure.snp.df, file='covered.4.1000.sure.snp.df.strand.nonspecific.exact.wilcox.downsampled.JvA.180807.rds')
