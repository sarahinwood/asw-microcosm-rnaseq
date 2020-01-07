library(tximport)
library(data.table)
library(DESeq2)

gene2tx <- fread("data/asw_mh_concat_transcriptome/ASW_Trinity.fasta.gene_trans_map", header = FALSE)
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

##create dds_concat_asw object and link to sample data  
dds_concat_asw <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##save dds_concat_asw object
saveRDS(dds_concat_asw, file = "output/concat_deseq2/dds_concat_asw.rds")

##create dds_concat_asw object for group analysis
dds_concat_group_asw <- copy(dds_concat_asw)

#######

##create groupings of weevil location and behavioural response to parasitoid
dds_concat_group_asw$group <- factor(paste(dds_concat_asw$Weevil_Location,dds_concat_asw$Behaviour,sep="_"))

##add group to design
design(dds_concat_group_asw) <- ~group
##run deseq2 and generate results
dds_concat_group_asw <- DESeq(dds_concat_group_asw)
##save dds_concat_asw_group
saveRDS(dds_concat_group_asw, file = "output/concat_deseq2/dds_concat_group_asw.rds")
