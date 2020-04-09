library("tximport")
library("data.table")
library("DESeq2")

gene2tx <- fread("data/asw_edited_transcript_ids/Trinity.fasta.gene_trans_map", header = FALSE)
tx2gene <- data.frame(gene2tx[, .(V2, V1)])

##Find all salmon quant files
quant_files <- list.files(path="output/asw_mh_concat_salmon/", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID - required for gene-level summarisation... 
##for methods that only provide transcript level estimaates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE)
##Import table describing samples
sample_data <- fread("data/sample_key_PCR_COI_counts.csv", header=TRUE)
setkey(sample_data, Sample_name)

##create dds object and link to sample data  
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)

dds_pc1 <- copy(dds)
##create sample groupings
dds_pc1$group <- factor(paste(dds$PC1_sign))

##add group to design
design(dds_pc1) <- ~group
##run deseq2 and generate results
dds_pc1 <- DESeq(dds_pc1)

resultsNames(dds_pc1)

pc1_res <- results(dds_pc1, contrast = c("group", "positive", "negative"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
pc1_ordered_res_group <- pc1_res[order(pc1_res$padj),]
##Make data table and write to output
pc1_ordered_res_group_table <- data.table(data.frame(pc1_ordered_res_group), keep.rownames = TRUE)
pc1_ordered_sig_res_group_table <- subset(pc1_ordered_res_group_table, padj < 0.05)
##merge with annots
asw_trinotate <- fread("data/asw_edited_transcript_ids/trinotate_longest_isoform.csv")
ordered_sig_degs_annots <- merge(pc1_ordered_sig_res_group_table, asw_trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE, all.y=FALSE)


plotCounts(dds_pc1, "ASW_TRINITY_DN2630_c0_g1", intgroup = c("group"), main="")


##fgsea analysis
trinotate_report <- asw_trinotate
gene_ids <- trinotate_report[!is.na(gene_ontology_Pfam), unique(`#gene_id`)]
full_res <- pc1_ordered_res_group_table

go_annot_list<-data.table(trinotate_report[,unique(unlist(strsplit(gene_ontology_Pfam, "`")))])
go_annot_table <- go_annot_list[,tstrsplit(V1, "^", fixed=TRUE)]
go_annot_table<-setnames(go_annot_table, old=c("V1", "V2", "V3"), new=c("pathway", "pathway_kind", "pathway_name"))

##function to extract GO terms from annotations in transcriptome (get all unique GO terms for each gene id) --> could look at other functional annot if I want to
EXTRACT_GO_TERMS <- function(x, trinotate_report){
  my_terms<-trinotate_report[`#gene_id`==x,unique(unlist(strsplit(gene_ontology_Pfam, "`")))]
  my_accessions<-unique(gsub("\\^.*", "", my_terms))
  my_accessions<-my_accessions[!is.na(my_accessions)]
  return(data.table(gene_id=x, accessions=my_accessions))
}

go_term_list <- lapply(gene_ids, EXTRACT_GO_TERMS, trinotate_report=trinotate_report)
go_term_table <- rbindlist(go_term_list)
term_to_gene <- go_term_table[,list(list(gene_id)), by=accessions]
pathways <- term_to_gene[,V1]
names(pathways) <- term_to_gene[,accessions]

##use stat column from deseq results to rank genes (can change if wanted)
setorder(full_res, stat)
ranks <- full_res[!is.na(stat), stat]
names(ranks) <- full_res[!is.na(stat), rn]

fgsea_res <- fgsea(pathways, ranks, nperm = 10000)
sorted_fgsea_res <- fgsea_res[order(fgsea_res$padj)]
##
sum(sorted_fgsea_res$padj<0.05)
## no sig enriched terms - seems unlikely/strange given 2000 DEGs??