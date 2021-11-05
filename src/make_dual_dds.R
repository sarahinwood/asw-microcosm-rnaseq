library(tximport)
library(tidyverse)
library(data.table)
library(DESeq2)

gene_trans_map <- 'data/asw-mh-combined-transcriptome/output/asw_mh_transcriptome/asw_mh_Trinity.fasta.gene_trans_map'
gene2tx <- fread(gene_trans_map, header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files
quant_files <- list.files(path="output/asw_mh_concat_salmon/", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimaates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_data <- fread("data/sample_table.csv", header=TRUE)
setkey(sample_data, sample_name)

##create dds object and link to sample data  
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##estimate size factors on whole dataset
dds <- estimateSizeFactors(dds)
##save whole dds object
saveRDS(dds, "output/deseq2/dual_dds.rds")

###############################################
##subset gene table into ASW, Mh & MhV genes ##
###############################################

##asw gene list
asw_tx <- subset(tx2gene, grepl("ASW_", V1))
asw_gene <- unique(asw_tx$V1)
##MhV gene list
mhv_genes <- fread("data/mh-rnaseq/output/blast/viral_genes/viral_genes_best_hits.csv")
mhv_genes$full_id <- paste("MH", mhv_genes$Trinity_ID, sep="_")
mhv_genes$full_id <- tstrsplit(mhv_genes$full_id, "_i", keep=c(1))
mhv_gene <- unique(mhv_genes$full_id)
##mh gene list
mh_tx <- subset(tx2gene, grepl("MH_", V1))
mh_gene_not_v <- setdiff(mh_tx$V1, mhv_gene)
mh_gene <- unique(mh_gene_not_v)

##subset asw dds
asw_dds <- dds[asw_gene,]
saveRDS(asw_dds, "output/deseq2/asw_dual/asw_dual_dds.rds")
##subset mh dds
mh_dds <- dds[mh_gene,]
saveRDS(mh_dds, "output/deseq2/mh_dual/mh_dual_dds.rds")
##subset MhV dds
mhv_dds <- dds[mhv_gene,]
saveRDS(mhv_dds, "output/deseq2/mhv_dual/mhv_dual_dds.rds")


