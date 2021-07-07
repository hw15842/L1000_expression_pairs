
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

# combine the cis and trans data 
trans_and_cis_combined <- rbind.fill(cis_data_trans_snps_only, sig_trans_eQTLs)

# add another new fdr column for all of them (not sure which one I should be using but being stringent and using this one but can go back and change)
trans_and_cis_combined$fdr_cis_trans_combined <- p.adjust(trans_and_cis_combined$Pvalue, "fdr")

# extract the significant ones based in the fdr pval from combined data
trans_and_cis_combined_sig <- subset(trans_and_cis_combined, trans_and_cis_combined$fdr_cis_trans_combined < 0.05)


### Split the data into dataframes based on the SNP

func_split <- function(SNP){

	snp_name <- trans_unique_eQTLs[SNP]
	print(snp_name)
	df <- subset(trans_and_cis_combined_sig, trans_and_cis_combined_sig$SNP == paste0(snp_name))
	return(df)
}


trans_and_cis_list_dfs <- lapply(1:length(trans_unique_eQTLs, func_split)



### create a list of the genes names against the other gene names within each dataframe(split by SNP already)


func_list <- function(dataframe){

	df <- trans_and_cis_list_dfs[[dataframe]]
	gene_list <- as.character(df$GeneSymbol)
	gene_pairs <- expand.grid(gene_list, gene_list)
	gene_pairs <- subset(gene_pairs, !gene_pairs$Var1 == gene_pairs$Var2)
	return(gene_pairs)

}


gene_pairs_list_dfs <- lapply(1:length(trans_and_cis_list_dfs), func_list)
gene_expression_pairs <- ldply(gene_pairs_list_dfs, data.frame)
gene_expression_pairs <- unique(gene_expression_pairs)

save(gene_expression_pairs, file="gene_expression_pairs.rdata")


















