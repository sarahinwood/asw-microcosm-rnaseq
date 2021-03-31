library("data.table")
library("DESeq2")
library("ggplot2")
library("RColorBrewer")

##finds genes not clearly DE jsut based off parasitism - stick with WT
##differences in gene expression based on ASW parasitism
##controlling for other factors

dds <- readRDS("output/deseq2/asw/asw_dds.rds")
dds_parasitism <- copy(dds)

dds_parasitism$Location <- factor(dds_parasitism$Weevil_Location)
dds_parasitism$Parasitism <- factor(dds_parasitism$Parasitism_status)
dds_parasitism$Sample_cleaned <- factor(dds_parasitism$Cleaned)

##Look for genes DE between parasitisms, controlling for parasitism status
design(dds_parasitism) <- ~Sample_cleaned+Location+Parasitism
dds_parasitism <- DESeq(dds_parasitism, test='LRT', reduced=~Sample_cleaned+Location)
saveRDS(dds_parasitism, "output/deseq2/asw/LRT_parasitism/dds_parasitism.rds")

dds_parasitism <- readRDS("output/deseq2/asw/LRT_parasitism/dds_parasitism.rds")

##extract results
dds_res <- results(dds_parasitism, alpha=0.1)

##results for all genes - for FGSEA analysis
full_res_dt <- data.table(data.frame(dds_res), keep.rownames=TRUE)
fwrite(full_res_dt, "output/deseq2/asw/LRT_parasitism/LRT_parasitism_para_full_res.csv")

##make list of sig genes
sig_degs <- subset(dds_res, padj<0.05)
sig_degs_names <- row.names(sig_degs)
fwrite(data.table(sig_degs_names), "output/deseq2/asw/LRT_parasitism/sig_degs_names.csv")

##sig degs table
ordered_sig_degs <- sig_degs[order(sig_degs$padj),]
ordered_sig_degs_table <- data.table(data.frame(ordered_sig_degs), keep.rownames = TRUE)
fwrite(ordered_sig_degs_table, "output/deseq2/asw/LRT_parasitism/sig_degs.csv")
##merge with annots
longest_trinotate_report <- fread("data/asw_longest_isoform_annots.csv")
ordered_sig_degs_annots <- merge(ordered_sig_degs_table, longest_trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE, all.y=FALSE)
fwrite(ordered_sig_degs_annots, "output/deseq2/asw/LRT_parasitism/sig_degs_annots.csv")

##plot counts for genes of interest, sub in name
plotCounts(dds_parasitism, "ASW_TRINITY_DN12172_c3_g1", intgroup = c("Parasitism", "parasitism"))

##can we do some sort of clustering based on expression pattern like in timecourse?
##do we really want to control for parasitism status though - we want differences between unpara ru and unpara dun


#########
##FGSEA##
#########
trinotate_report <- fread("data/asw_edited_transcript_ids/trinotate_annotation_report.txt", na.strings=".")
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
setorder(full_res_dt, stat)
ranks <- full_res_dt[!is.na(stat), stat]
names(ranks) <- full_res_dt[!is.na(stat), rn]

fgsea_res <- fgsea(pathways, ranks, nperm = 10000)
sorted_fgsea_res <- fgsea_res[order(fgsea_res$padj)]
##44 enriched GO terms
sum(sorted_fgsea_res$padj<0.05)
sig_fgsea_res <- subset(sorted_fgsea_res, padj < 0.05)
annot_sig_fgsea <- merge(sig_fgsea_res, go_annot_table, by.x="pathway", by.y="pathway", all.x=TRUE)
fwrite(annot_sig_fgsea, "output/deseq2/asw/LRT_parasitism/sig_annot_fgsea_pfam.csv")

##split into 3 tables --> biological process, cellular component and molecular function
bp_res <- annot_sig_fgsea[annot_sig_fgsea$pathway_kind=="biological_process"]
cc_res <- annot_sig_fgsea[annot_sig_fgsea$pathway_kind=="cellular_component"]
mf_res <- annot_sig_fgsea[annot_sig_fgsea$pathway_kind=="molecular_function"]

##plot normalised enrichment for GO terms where padj<0.1 (but indicate if padj<0.05) - can change to only bp, cc or mf
ggplot(mf_res, aes(reorder(pathway_name, NES), NES)) +
  geom_text(aes(label=round(padj, digits=3)), vjust=0, hjust=0) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Molecular Function GO Pathway", y="FGSEA Normalized Enrichment Score") + 
  theme_minimal()+
  theme(axis.text.y = element_text(size=15))


