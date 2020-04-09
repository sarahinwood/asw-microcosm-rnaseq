library(data.table)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(EnhancedVolcano)
library(fgsea)

dds_location_attacked <- readRDS("output/deseq2/asw/WT_location_attacked/dds_location_attacked.rds")

resultsNames(dds_location_attacked)

##Make table of results for exposed vs control heads
res_group <- results(dds_location_attacked, contrast = c("group", "Ruakura_Y", "Dunedin_Y"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
##write tables
fwrite(ordered_res_group_table, "output/deseq2/asw/WT_location_attacked/location_para/full_res.csv")
fwrite(ordered_sig_res_group_table, "output/deseq2/asw/WT_location_attacked/location_para/sig_degs.csv", col.names = TRUE, row.names = FALSE)

###merge with annots
longest_trinotate_report <- fread("data/asw_edited_transcript_ids/trinotate_longest_isoform.csv")
sig_degs_annots <- merge(ordered_sig_res_group_table, longest_trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_degs_annots, "output/deseq2/asw/WT_location_attacked/location_para/sig_degs_annots.csv")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", pointSize = 3, pCutoff=0.05)

plotCounts(dds_location_attacked, "ASW_TRINITY_DN15587_c0_g1", intgroup = c("group"))

##################
##FGSEA analysis##
##################
trinotate_report <- fread("data/asw_edited_transcript_ids/trinotate_annotation_report.txt", na.strings=".")
full_res <- ordered_res_group_table
gene_ids <- trinotate_report[!is.na(gene_ontology_Pfam), unique(`#gene_id`)]

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
fwrite(sorted_fgsea_res, "output/deseq2/asw/WT_location_attacked/location_para/fgsea_pfam_res.csv")

sig_fgsea_res <- subset(sorted_fgsea_res, padj < 0.1)
annot_sig_fgsea <- merge(sig_fgsea_res, go_annot_table, by.x="pathway", by.y="pathway", all.x=TRUE)
fwrite(annot_sig_fgsea, "output/deseq2/asw/WT_location_attacked/location_para/sig_annot_fgsea_pfam.csv")
