library(tximport)
library(data.table)
library(DESeq2)

gene2tx <- fread("data/asw_mh_concat_transcriptome/MH_Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files from concat asw_mh mapping
quant_files <- list.files(path="output/asw_mh_concat_salmon/", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_data <- fread("data/sample_key.csv", header=TRUE)
setkey(sample_data, Sample_name)

##create dds_concat_mh object and link to sample data  
dds_concat_mh <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##save dds_concat_mh object
saveRDS(dds_concat_mh, file = "output/concat_deseq2/dds_concat_mh.rds")

##create dds_concat_mh object for group analysis
dds_concat_group_mh <- copy(dds_concat_mh)

#######

##create groupings of weevil location and behavioural response to parasitoid
dds_concat_group_mh$group <- factor(paste(dds_concat_mh$Behaviour,sep="_"))

##add group to design
design(dds_concat_group_mh) <- ~group
##run deseq2 and generate results
dds_concat_group_mh <- DESeq(dds_concat_group_mh)
##save dds_concat_group_mh
saveRDS(dds_concat_group_mh, file = "output/concat_deseq2/dds_concat_group_mh.rds")
