#!/usr/bin/env Rscript

#this script is used to, per chromosome:
#DOWNSAMPLE the data
#REFORMAT the data to facilitate downstream analysis
#NORMALIZE the data and combine the different libraries and replicates



# retrieve argument from commandline. User has to supply a single chromosome name as arg.
args = commandArgs(trailingOnly=TRUE)
print(args[1])

# Change to directory which will be used for all output files
setwd("~/mydata/git/SuRE_shared/Joris/SNP_project/all_snp_analysis_wilcox/")

# The vector with all chromosome names
chr.vector <- c(paste('chr', seq(1:22), sep=''),'chrX')


#install.packages("splitstackshape")
library(splitstackshape)

# The directories containing the SuRE-count tables per library
# The names of the files for chr1 are given. In below code "chr1" is replaced
# by the chromosome name supplied by the user.
path.vector <- c(
  "/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE42/SuRE-pipelineOutput/SuRE-counts_chr1.inf.txt.gz",
  "/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP180619_SuRE42_2/pipelineOutput/SuRE-counts_chr1.txt.gz",
  "/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE43/SuRE43-pipelineOutput/SuRE-counts_chr1.inf.txt.gz",
  "/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP180619_SuRE43_2/pipelineOutput/SuRE-counts_chr1.txt.gz",
   "/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE44/SuRE44-pipelineOutput/SuRE-counts_chr1.inf.txt.gz",
"/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP180619_SuRE44_2/pipelineOutput/SuRE-counts_chr1.txt.gz",
"/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE45_1/SuRE45_1-pipelineOutput/SuRE-counts_chr1.inf.txt.gz",
"/DATA/usr/ludo/projects/LP140430_SureSeq_JvArensbergen/analyses/LP170216_SuRE45_2/SuRE45_2-pipelineOutput/SuRE-counts_chr1.inf.txt.gz")
                 
#Get the previously determined fractions to which samples should be downsampled
k562.downsampling.vector.SuRE42_1 <- readRDS('k562.downsampling.vector.SuRE42_1.180803.rds')
HepG2.downsampling.vector.SuRE42_1 <- readRDS('HepG2.downsampling.vector.SuRE42_1.180803.rds')
HepG2.downsampling.vector.SuRE43_1 <- readRDS('HepG2.downsampling.vector.SuRE43_1.180803.rds')


# convert the user argument (ie chromosome name) to integer
k <- as.numeric(args[1])

## START OF DOWNSAMPLING ##
paste('now started analysis of ',chr.vector[k],' (counting up)', sep='')

sure42_1.snp <- fread(paste('zcat',gsub('chr1',chr.vector[k],path.vector[1])), sep='\t', header=TRUE,stringsAsFactors=FALSE, data.table=FALSE) 
sure42_2.snp <- fread(paste('zcat',gsub('chr1',chr.vector[k],path.vector[2])), sep='\t', header=TRUE,stringsAsFactors=FALSE, data.table=FALSE)
sure43_1.snp <- fread(paste('zcat',gsub('chr1',chr.vector[k],path.vector[3])), sep='\t', header=TRUE,stringsAsFactors=FALSE, data.table=FALSE) 
sure43_2.snp <- fread(paste('zcat',gsub('chr1',chr.vector[k],path.vector[4])), sep='\t', header=TRUE,stringsAsFactors=FALSE, data.table=FALSE) 
sure44_1.snp <- fread(paste('zcat',gsub('chr1',chr.vector[k],path.vector[5])), sep='\t', header=TRUE,stringsAsFactors=FALSE, data.table=FALSE) 
sure44_2.snp <- fread(paste('zcat',gsub('chr1',chr.vector[k],path.vector[6])), sep='\t', header=TRUE,stringsAsFactors=FALSE, data.table=FALSE) 
sure45_1.snp <- fread(paste('zcat',gsub('chr1',chr.vector[k],path.vector[7])), sep='\t', header=TRUE,stringsAsFactors=FALSE, data.table=FALSE) 
sure45_2.snp <- fread(paste('zcat',gsub('chr1',chr.vector[k],path.vector[8])), sep='\t', header=TRUE,stringsAsFactors=FALSE, data.table=FALSE) 
 

##need to switch columns around for (SuRE42_2, 43_2, 44_2). These were separately processed and column orders were not the same as for previous samples
colnames(sure42_1.snp)
colnames(sure42_2.snp)

sure42_2.snp <- sure42_2.snp[,c(1:10,13:15,11,12,16:19)]
sure43_2.snp <- sure43_2.snp[,c(1:10,13:15,11,12,16:19)]
sure44_2.snp <- sure44_2.snp[,c(1:10,13:15,11,12,16:19)]

# Start with downsampling SuRE42_1 and SuRE43_1, as these samples are sequenced much deeper than the others
##for sampling, give each row for SuRE42_1 and SuRE43_1 a row identifier
sure42_1.snp$number <- seq(1:nrow(sure42_1.snp))
sure43_1.snp$number <- seq(1:nrow(sure43_1.snp))

# downsample all 5 cDNA samples (3x K562, 2x HepG2) for SuRE42_1
for(i in 1:5){
  snpcounts <- rep(sure42_1.snp$number,sure42_1.snp[,i+10])##simply make a big vector of the counts with the rowidentifiers from which you can then sample
  head(snpcounts)
  
  snpcounts.sample <- sample(snpcounts, size=round(length(snpcounts)*c(k562.downsampling.vector.SuRE42_1,HepG2.downsampling.vector.SuRE42_1)[i],digits=0))
  
  sure42_1.snp[,i+10] <- 0
  sure42_1.snp[match(names(table(snpcounts.sample)),sure42_1.snp$number),i+10] <- table(snpcounts.sample)
}

# downsample all 2 HepG2 cDNA samples for SuRE43_1
for(i in 1:2){
  snpcounts <- rep(sure43_1.snp$number,sure43_1.snp[,13+i])
  head(snpcounts)
  snpcounts.sample <- sample(snpcounts, size=round(length(snpcounts)*HepG2.downsampling.vector.SuRE43_1[i],digits=0))
  
  sure43_1.snp[,i+13] <- 0
  sure43_1.snp[match(names(table(snpcounts.sample)),sure43_1.snp$number),i+13] <- table(snpcounts.sample)

}

# several checks to make sure sampling was correct:
head(sure42_1.snp)
head(sure42_2.snp)

#do they now have similar sums?
colSums(sure42_1.snp[,11:15])
colSums(sure43_1.snp[,11:15])
colSums(sure42_2.snp[,11:15]) 


##remove the last column containing the row identifying used for subsampling
sure42_1.snp <- sure42_1.snp[,1:18]
sure43_1.snp <- sure43_1.snp[,1:18]

# Write the downsampled data to file:
write.table(sure42_1.snp, file=paste0('~/mydata/git_data/Downsampled_SuRE_counts_180802/SuRE42_1_downsampled_',chr.vector[k],'.txt'), quote=FALSE, col.names=NA)
write.table(sure43_1.snp, file=paste0('~/mydata/git_data/Downsampled_SuRE_counts_180802/SuRE43_1_downsampled_',chr.vector[k],'.txt'), quote=FALSE, col.names=NA)
## End of DOWNSAMPLING data ##


## start REFORMAT the data to simplify subsequent processing ##
# current data may contain multiple SNP annotations per fragment.
# Here these cases are separated to have separated rows for each individual SNP

# collect all data in a list
sure.list <- list(sure42_1.snp,sure42_2.snp,sure43_1.snp,sure43_2.snp,sure44_1.snp,sure44_2.snp,sure45_1.snp,sure45_2.snp)

# discard all rows not containing any SNP:
sure.list.reduced <- list()
for(i in 1: length(sure.list)){
  sure.list.reduced [[i]] <- sure.list[[i]][which(!sure.list[[i]]$SNPidx==""),]
}
print('finished reducing to SNP containing')

# Columns with SNP annotation contain comma separated values in case multiple SNP's overlap a certain fragment
# Here these values are split into separate rows, duplicating all other columns
# cSplit(..) from package 'splitstackshape' is used to split the column with
# SNP absolute positions. The values are split in vertical (ie "long")
# direction, resulting in additional rows for every split.

sure.split.list <- list()
for(i in 1:length(sure.list.reduced)){
  print(i)
  sure.split.list[[i]] <- cSplit(sure.list.reduced[[i]][,c(1:4,8,10:15)], c("SNPabspos"), sep = ",", direction = "long") 
}
# check
sure.split.list[[1]][1:50,]

# add the columns for ref/alternative annotation and base
for(i in 1: length(sure.split.list)){
  print(i)
  sure.split.list[[i]][,'SNPvarInf'] <- unlist(strsplit(sure.list.reduced[[i]][,'SNPvarInf'],","))
  sure.split.list[[i]][,'SNP_ID'] <- unlist(strsplit(sure.list.reduced[[i]][,'SNP_ID'],","))
}
print('finished splitting')
# inspect the splitting was succesful
head(sure.split.list[[1]])
dim(sure.split.list[[1]])
class(sure.split.list[[1]])

# coerce the data.table's into dataframes
for(i in 1: length(sure.split.list)){
  print(i)
  sure.split.list[[i]] <- as.data.frame(sure.split.list[[i]])
}

#there are duplicated rows where a SNP was read 2 times (i.e. from both ends), because SNPvar column then contains the same SNP twice. 
#get rid of those duplciated rows. 

for(i in 1: length(sure.split.list)){
  print(i)
  sure.split.list[[i]] <- sure.split.list[[i]][!duplicated(sure.split.list[[i]]),]
}

rm(sure.list.reduced)
gc()
print('finished removing duplicates')
rm(sure.list,sure42_1.snp,sure42_2.snp,sure43_1.snp,sure43_2.snp,sure44_1.snp,sure44_2.snp,sure45_1.snp,sure45_2.snp)
# done REFORMAT data 


## start NORMALIZE counts ## 
for(i in 1: length(sure.split.list)){
  print(i)
  # assign column names
  colnames(sure.split.list[[i]])[7:11] <- c('K562.cDNA1','K562.cDNA2','K562.cDNA3','HEPG2.cDNA1','HEPG2.cDNA2')


  # express cDNA as read per billion, rounded to whole integers
  sure.split.list[[i]]$norm.K562.cDNA1 <- round((1e9/sum(sure.split.list[[i]]$K562.cDNA1)) * sure.split.list[[i]]$K562.cDNA1,digits=0)
  sure.split.list[[i]]$norm.K562.cDNA2 <- round((1e9/sum(sure.split.list[[i]]$K562.cDNA2)) * sure.split.list[[i]]$K562.cDNA2,digits=0)
  sure.split.list[[i]]$norm.K562.cDNA3 <- round((1e9/sum(sure.split.list[[i]]$K562.cDNA3)) * sure.split.list[[i]]$K562.cDNA3,digits=0)
  sure.split.list[[i]]$norm.HEPG2.cDNA1 <- round((1e9/sum(sure.split.list[[i]]$HEPG2.cDNA1)) * sure.split.list[[i]]$HEPG2.cDNA1,digits=0)
  sure.split.list[[i]]$norm.HEPG2.cDNA2 <- round((1e9/sum(sure.split.list[[i]]$HEPG2.cDNA2)) * sure.split.list[[i]]$HEPG2.cDNA2,digits=0)
  # normalize iPCR
  sure.split.list[[i]]$norm.iPCR <- round((1e9/sum(sure.split.list[[i]]$iPCR)) * sure.split.list[[i]]$iPCR,digits=0)
}
# inspect results:
sure.split.list[[1]][1:10,]

# add a column to indicate library
sure.split.list[[1]]$library <- 'SuRE42_1'
sure.split.list[[2]]$library <- 'SuRE42_2'
sure.split.list[[3]]$library <- 'SuRE43_1'
sure.split.list[[4]]$library <- 'SuRE43_2'
sure.split.list[[5]]$library <- 'SuRE44_1'
sure.split.list[[6]]$library <- 'SuRE44_2'
sure.split.list[[7]]$library <- 'SuRE45_1'
sure.split.list[[8]]$library <- 'SuRE45_2'

# combine the data over all libraries into a sine dataframe
sure.snp.combined.df <- rbind(sure.split.list[[1]],sure.split.list[[2]],sure.split.list[[3]],sure.split.list[[4]],sure.split.list[[5]],sure.split.list[[6]],sure.split.list[[7]],sure.split.list[[8]])

# Combine the cDNA counts of replicate transfections, and divide by number of replicates 
sure.snp.combined.df$sum.norm.K562.cDNA <- (sure.snp.combined.df$norm.K562.cDNA1+sure.snp.combined.df$norm.K562.cDNA2+sure.snp.combined.df$norm.K562.cDNA3)/3
sure.snp.combined.df$sum.norm.HEPG2.cDNA <- (sure.snp.combined.df$norm.HEPG2.cDNA1+sure.snp.combined.df$norm.HEPG2.cDNA2)/2

# Normalize cDNA counts relative to iPCR counts
sure.snp.combined.df$ipcr.norm.sum.K562.cDNA <- sure.snp.combined.df$sum.norm.K562.cDNA/sure.snp.combined.df$norm.iPCR
sure.snp.combined.df$ipcr.norm.sum.HEPG2.cDNA <- sure.snp.combined.df$sum.norm.HEPG2.cDNA/sure.snp.combined.df$norm.iPCR

# inspect:
head(sure.snp.combined.df)

##split the plus and minus strand
sure.snp.combined.plus.df <- sure.snp.combined.df[sure.snp.combined.df$strand=='+',]
sure.snp.combined.minus.df <- sure.snp.combined.df[sure.snp.combined.df$strand=='-',]

# save the data
saveRDS(sure.snp.combined.plus.df, file=paste('~/mydata/git_data/Downsampled_SuRE_counts_180802/combined.sure.snp.files.180802/sure.snp.combined.plus.',chr.vector[k],'.',Sys.Date(),'.rds', sep=''))
saveRDS(sure.snp.combined.minus.df, file=paste('~/mydata/git_data/Downsampled_SuRE_counts_180802/combined.sure.snp.files.180802/sure.snp.combined.minus.',chr.vector[k],'.',Sys.Date(),'.rds', sep=''))

paste('now finished analysis of ',chr.vector[k],' (counting up)', sep='')

