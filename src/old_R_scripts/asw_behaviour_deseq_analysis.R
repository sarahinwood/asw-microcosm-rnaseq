library("tximport")
library("data.table")
library("DESeq2")

gene2tx <- fread("data/asw_transcriptome/Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files from salmon filtering because res. from filtering with STAR and lost bro that way
quant_files <- list.files(path="output/asw_salmon/", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_data <- fread("data/sample_key.csv", header=TRUE)
setkey(sample_data, Sample_name)

##create dds object and link to sample data  
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
##save dds object
saveRDS(dds, file = "output/deseq2/dds.rds")

##create dds object for group analysis
dds_group <- copy(dds)

#######

##create groupings of weevil location and behavioural response to parasitoid
dds_group$group <- factor(paste(dds$Weevil_Location,dds$Behaviour,sep="_"))

##add group to design
design(dds_group) <- ~group
##run deseq2 and generate results
dds_group <- DESeq(dds_group)
##save dds_group
saveRDS(dds_group, file = "output/deseq2/dds_group.rds")