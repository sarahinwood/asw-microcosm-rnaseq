library(tximport)
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

##subset gene table into ASW and Mh genes
asw_tx <- subset(tx2gene, grepl("ASW_", V1))
asw_gene <- unique(asw_tx$V1)
mh_tx <- subset(tx2gene, grepl("MH_", V1))
mh_gene <- unique(mh_tx$V1)

##subset dds
asw_dds <- dds[asw_gene,]
saveRDS(asw_dds, "output/deseq2/asw_dual/asw_dual_dds.rds")

mh_dds <- dds[mh_gene,]
saveRDS(mh_dds, "output/deseq2/mh_dual/mh_dual_dds.rds")

##save whole dds object
saveRDS(dds, "output/deseq2/dual_dds.rds")
