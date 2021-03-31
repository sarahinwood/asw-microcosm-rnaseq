library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(fgsea)

##Pairwise between para and para-undetected

dds <- readRDS("output/deseq2/mh/mh_dds.rds")

##create dds object for parasitism analysis
dds_parasitism <- copy(dds)
dds_parasitism$group <- factor(paste(dds$Parasitism_status))

##add group to design
design(dds_parasitism) <- ~group
##run deseq2 and generate results
dds_parasitism <- DESeq(dds_parasitism)
##save dds_group
saveRDS(dds_parasitism, file = "output/deseq2/mh/WT_parasitism/mh_dds_parasitism.rds")

dds_parasitism <- readRDS("output/deseq2/mh/WT_parasitism/mh_dds_parasitism.rds")

resultsNames(dds_parasitism)

##Make table of results for exposed vs control heads
res_group <- results(dds_parasitism, contrast = c("group", "parasitised", "undetected"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
##write tables
fwrite(ordered_res_group_table, "output/deseq2/mh/WT_parasitism/full_res.csv")
fwrite(ordered_sig_res_group_table, "output/deseq2/mh/WT_parasitism/sig_degs.csv", col.names = TRUE, row.names = FALSE)

mh_viral_ids <- fread("data/mh_viral_transcript_ids.txt")
mh_viral_ids$fixed_id <- paste("MH", mh_viral_ids$`#gene_id`, sep="_")
viral_transcript_results <- subset(ordered_res_group_table, rn %in% mh_viral_ids$fixed_id)
viral_transcript_results$sig <- (viral_transcript_results$padj<0.05)
mh_viral_expression <- fread("data/mh_viral_expression_matrix.csv")
mh_viral_expression$rn <- paste("MH", mh_viral_expression$`#gene_id`, sep="_")
mh_viral_transcripts_res <- merge(mh_viral_expression, viral_transcript_results, by="rn")
mh_viral_transcripts_res_sig <- subset(mh_viral_transcripts_res, sig=="TRUE")
fwrite(mh_viral_transcripts_res_sig, "output/deseq2/mh/WT_parasitism/viral_degs_mh_expression.csv")

###merge with annots
longest_trinotate_report <- fread("data/mh_edited_transcript_ids/trinotate_longest_isoform.csv")
sig_degs_annots <- merge(ordered_sig_res_group_table, longest_trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_degs_annots, "output/deseq2/mh/WT_parasitism/sig_degs_annots.csv")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", pointSize = 3, pCutoff=0.05)

##can add in parasitism to check DE isn't a result of parasitism
plotCounts(dds_parasitism, "MH_TRINITY_DN1053_c1_g1", intgroup = c("group"))

plot_gene <- plotCounts(dds_parasitism, "MH_TRINITY_DN11733_c0_g1", 
                        intgroup = c("group"), returnData = TRUE)
ggplot(plot_gene,
       aes(x = group, y = count)) + 
  geom_point() + geom_smooth(se = FALSE, method = "loess")+
  scale_y_log10() + ylab("Normalised Count")+
  xlab("Parasitism Status")+
  ggtitle("MH_TRINITY_DN11733_c0_g1 - KilA")

#########
##FGSEA##
#########
trinotate_report <- fread("data/mh_edited_transcript_ids/trinotate_annotation_report.txt", na.strings = ".")
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
setorder(ordered_res_group_table, stat)
ranks <- ordered_res_group_table[!is.na(stat), stat]
names(ranks) <- ordered_res_group_table[!is.na(stat), rn]

fgsea_res <- fgsea(pathways, ranks, nperm = 10000)
sorted_fgsea_res <- fgsea_res[order(fgsea_res$padj)]
sum(sorted_fgsea_res$padj<0.05)

sig_fgsea_res <- subset(sorted_fgsea_res, padj < 0.1)
annot_sig_fgsea <- merge(sig_fgsea_res, go_annot_table, by.x="pathway", by.y="pathway", all.x=TRUE)
fwrite(annot_sig_fgsea, "output/deseq2/mh/WT_parasitism/fgsea_annot_sig_terms.csv")

