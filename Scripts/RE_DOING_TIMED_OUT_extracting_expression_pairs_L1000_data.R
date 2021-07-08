

################################################
### RE_DOING_TIMED_OUT_extracting_expression_pairs_L1000_data.R ###
################################################

## module add languages/r/4.0.3 
# has to be in R version 4

args  <-  commandArgs(trailingOnly=TRUE)
data_location <- toString(args[1])


setwd(paste0(data_location))

library("cmapR")
library("data.table")
library("dplyr")
library("biomaRt")

load("gene_expression_pairs_l1000_list_of_dataframes.rdata")

gene_expression_pairs_l1000_results = do.call(rbind, datalist)


colnames <- c("prot_1", "prot_2", "sig_info","Z") 
datalist <- lapply(datalist, setNames, colnames)

colnames(gene_expression_pairs_l1000_results) <- c("gene_1", "gene_2", "sig_info","Z")

# save(datalist, file="gene_expression_pairs_l1000_list_of_dataframes.rdata")

save(gene_expression_pairs_l1000_results, file="gene_expression_pairs_l1000_results.rdata")

# datalist gives us a list of dataframes - one dataframe for each cpg snp pair - 3,390 dataframes
# prot_pairs_l1000_results gives us dataframe of all snp-cpg pairs and their associated sig_id and the z-score for that sig_id
# sig_id = (A CMap unique identification number assigned to each signature generated from L1000 data)


# Z-score   For a data point, the number of standard deviations that point is above or below the population mean is called its Z-score. 
# In the L1000 data processing pipeline, we compute a robust z-score for each gene in each sample. 
# The reference population used to compute the median and MAD is the expression of the given gene in every other well on the plate. 
# These z-score values correspond to level 4 data



### startby taking an average of all the different perturbations for each protein pair

ave_zscore_func <- function(gene_expression_pairs_dataset){

	df <- datalist[[gene_expression_pairs_dataset]]
	z_ave <- mean(df$Z)
	df <- df[1, 1:3]
	df <- as.data.frame(cbind(df, z_ave))
	names(df) <- c("prot_1", "prot_2", "sig_info", "Z_ave")
	return(df)
}

gene_expression_pairs_l1000_ave_zscores <- lapply(1:length(datalist), ave_zscore_func)
gene_expression_pairs_l1000_ave_zscores <- ldply(gene_expression_pairs_l1000_ave_zscores, data.table)

save(gene_expression_pairs_l1000_ave_zscores, file="gene_expression_pairs_l1000_ave_zscores.rdata")



