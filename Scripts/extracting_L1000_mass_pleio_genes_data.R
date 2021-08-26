
################################################
### extracting_L1000_mass_pleio_genes_data.R ###
################################################

## module add languages/r/4.0.3 
# has to be in R version 4

args  <-  commandArgs(trailingOnly=TRUE)
data_location <- toString(args[1])

# data_location <- "/mnt/storage/scratch/hw15842/repo/L1000_expression_pairs/Data/"

setwd(paste0(data_location))

library("cmapR")
library("data.table")
library("dplyr")
library("biomaRt")
library("plyr")


## First part is exactly the same as Tom has done - only looking at certain perturbations as well

## read in the moderated Z scores file ##
# rid = either a vector of character or integer row indices or a path to a grp file containing character row indices. Only these indices will be parsed from the file.
# not sure why rid = 1 here .....
ds <- parse.gctx(paste0(data_location, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"), rid=1)

# mat - the data matrix
# rdesc - a data.frame of row annotations, with one row per matrix row
# cdesc - a data.frame of column annotations, with one row per matrix column
# rid - a character vector of unique row identifiers
# cid - a character vector of unique column identifiers
# src - a character string indicating the source (usually a file path) of the data


## In the ds file, extract the cids section
cids <- as.data.frame(ds@cid)

# call the column names of the cid dataframe "sig_name"
colnames(cids) <- c('sig_name')

# add a column to the cids dataframe called "index" which is just the row names (1 to 473647)
cids$index <- row.names(cids)

# read in the sig_info.txt file 
# "Metadata for each signature in the Level 5 matrix (metadata for the columns in the Level 5 data matrix)""
a <- fread(paste0(data_location, "GSE92742_Broad_LINCS_sig_info.txt"))

## filter the sig_infor.txt file based on the perturbagen type 
# here tom has kept:
# 	- trt_sh.cgs' = Consensus signature from shRNAs targeting the same gene
#	- 'trt_oe'  = cDNA for overexpression of wild-type gene
b <- filter(a, pert_type=='trt_sh.cgs'|pert_type=='trt_oe')

# subset the cid dataframe based on the sig_name column IDs that are in the filtered perturbagens file 
# takes it down to 58,925 rows (from 473,647 rows)
c <- cids[which(cids$sig_name %in% b$sig_id),]

# create a list of the index (row names) from datafrme c - i.e. the row names from the subsetsed pertrubagens list 
index <- as.numeric(c$index)

# select just the "sig_id" and "pert_iname" columns from the b dataset (the filtered pertrubagens dataset)
d <- b[,c('sig_id','pert_iname')]

# merge together d (sig_ id and pert_iname columns of the subsetsed pertrubagens dataset) and c (the )
d1 <- merge(d,c,by.x='sig_id',by.y='sig_name',all.x=T)
d2 <- d1[order(d1$index),]





##load data matrix of Z scores - re load the data but this time with just the "index" columns - the ones we have decided we want after filtering on the perturbations 
# again takes it from 473,647 data points to 58,925 columns
ds <- parse.gctx(paste0(data_location, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"), cid=index)



## map entrez genes to ENSG - think this is the more useful one according to Tom (ensemble # out in his script)
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
load(paste0(data_location, "ensembl_database_downloaded_from_biomart.R")) # loading it rather than trying to extract it from site as doesnt always work from site - however if running at a later data may need to see if anything has changed. 
f <- getBM(attributes=c('hgnc_symbol','entrezgene_id'), mart=ensembl) # The entrezgene attribute was changed to entrezgene_id when Ensembl 97 was released

# in the new ds that just has the purturbagens Tom wanted, extract the row names from the matrix 
entrez <- as.data.frame(row.names(ds@mat))
colnames(entrez) <- c("entrez_id")

# add an index column which is the row names of the entrez data frame
entrez$index <- row.names(entrez)

# merge the entrez ids from the ds and the ones we got from biomart, keeping all the ones from the ds - gives us the entrez ID, the index and the hgnc_symbol
merge_entrez <- merge(entrez,f,by.x="entrez_id","entrezgene_id",all.x=T)
# keep just the unique entrez_ids
uniq_entrez <- merge_entrez[!duplicated(merge_entrez$entrez_id),]
# order them numerically based on the index value (which is the row names of the entrz dataframe made from the row names of the ds matrix)
sort_entrez <- uniq_entrez[order(as.numeric(uniq_entrez$index)),]


# Load in the mass pleio genes 
# made in the Rmarkdown doc "L1000_analysis_using_expression_pairs"

load("mass_pleio_gene_list.rdata")

gene_list <- mass_pleio_gene_list$gene_1

func1 <- function(gene){

cols <- which(d2$pert_iname==as.character(paste0(gene)))

x_new <- colMeans(abs(ds@mat[,cols]))
dat_new <- tibble(id=names(x_new), z=x_new)
dat_new <- tidyr::separate(dat_new, id, sep=":", into=c("exp", "gene", "num"))
dat_new$gene_1 <- paste0(gene)
return(dat_new)

}


mass_pleio_genes_L1000_z_scores <- lapply(gene_list, func1)

save(mass_pleio_genes_L1000_z_scores, file="mass_pleio_genes_L1000_z_scores.rdata")


func2 <- function(dataset){

	df <- mass_pleio_genes_L1000_z_scores[[dataset]]
	mean_Z <- mean(df$z)
	dat <- data.frame(cbind(df[1,5],nrow(df),mean_Z))
	names(dat)[2] <- c("num_cols_gene")
	return(dat)
}


mass_pleio_genes_L1000_z_scores_summary <- lapply(1:length(mass_pleio_genes_L1000_z_scores), func2)
mass_pleio_genes_L1000_z_scores_summary <- bind_rows(mass_pleio_genes_L1000_z_scores_summary)

save(mass_pleio_genes_L1000_z_scores_summary, file="mass_pleio_genes_L1000_z_scores_summary.rdata")


