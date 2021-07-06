## extract_expression_pairs_eQTLgen.R ##


args  <-  commandArgs(trailingOnly=TRUE)
data_location <- toString(args[1])

# data_location <- "/mnt/storage/scratch/hw15842/repo/L1000_expression_pairs/Data/"

library("data.table")
library("dplyr")
library("plyr")

setwd(paste0(data_location))

sig_trans_eQTLs_filename <- gzfile('2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz','rt')  
sig_trans_eQTLs <- read.table(sig_trans_eQTLs_filename,header=T)

## ran the below once and then saved, took ages to load in so just subsetted the data we needed already and then can load the smaller file back in 

# full_cis_eQTLs_filename <- gzfile('2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz','rt')  
# full_cis_eQTLs <- read.table(full_cis_eQTLs_filename,header=T)
# trans_unique_eQTLs <- unique(sig_trans_eQTLs$SNP)
# cis_data_trans_snps_only <- subset(full_cis_eQTLs, full_cis_eQTLs$SNP %in% trans_unique_eQTLs)
# save(cis_data_trans_snps_only, file="cis_data_trans_snps_only.rdata")

load("cis_data_trans_snps_only.rdata")

trans_unique_eQTLs <- unique(sig_trans_eQTLs$SNP)


# create a new FDR column

cis_data_trans_snps_only$fdr_new <- p.adjust(cis_data_trans_snps_only$Pvalue, "fdr")


