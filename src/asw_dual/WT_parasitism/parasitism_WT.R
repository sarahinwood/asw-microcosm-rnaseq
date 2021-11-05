library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(viridis)
library(pheatmap)
library(tidyverse)
library(dplyr)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
##factors and design
asw_dds$para <- factor(paste(asw_dds$Parasitism_status))
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$pc1_sign <- factor(paste(asw_dds$PC1_sign))
##relevel parasitism
asw_dds$para <- relevel(asw_dds$para, ref="undetected")

asw_dds_para <- copy(asw_dds)
design(asw_dds_para) <- ~pc1_sign+location+para
##run deseq2
asw_dds_para <- DESeq(asw_dds_para)
saveRDS(asw_dds_para, "output/deseq2/asw_dual/WT_parasitism/parasitism_WT.rds")

asw_dds_para <- readRDS("output/deseq2/asw_dual/WT_parasitism/parasitism_WT.rds")
##results
res_group <- results(asw_dds_para, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/asw_dual/WT_parasitism/sig_degs.csv")
fwrite(ordered_res_group_table, "output/deseq2/asw_dual/WT_parasitism/res_group.csv")

plotCounts(asw_dds_para, "ASW_TRINITY_DN6332_c1_g1", intgroup="para")
##volcano plot
EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.1, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

##merge with trinotate annots
trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw_dual/WT_parasitism/sig_annots.csv")

##blastx for unann
unann_deg_blast <- fread("output/deseq2/asw_dual/WT_parasitism/blastx.outfmt6")
setnames(unann_deg_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("transcript_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
unann_annots <- unann_deg_blast[,c(1,11,13)]
unann_annots$transcript_id <- tstrsplit(unann_annots$transcript_id, "_i", keep=c(1))

##merge unann blast with all other annots
sig_all_annots <- merge(sig_annots, unann_annots, by.x="rn", by.y="transcript_id", all.x=TRUE)
fwrite(sig_all_annots, "output/deseq2/asw_dual/WT_parasitism/sig_all_annots.csv")

#############
## heatmap ##
#############
all_sig_annots <- fread("output/deseq2/asw_dual/WT_parasitism/sig_blast_annots.csv")
annotated_tbx <- subset(all_sig_annots, sprot_Top_BLASTX_hit!="")
annotated_nrbx <- subset(all_sig_annots, annotation!="")
all_annot <- full_join(annotated_tbx, annotated_nrbx)
##vst transform
asw_vst <- varianceStabilizingTransformation(asw_dds_para, blind=TRUE)
asw_vst_assay_dt <- data.table(assay(asw_vst), keep.rownames=TRUE)
##subset for DEGs
asw_vst_degs <- subset(asw_vst_assay_dt, rn %in% all_annot$rn)
asw_vst_degs$rn <- tstrsplit(asw_vst_degs$rn, "ASW_", keep=c(2))
##turn first row back to row name
asw_vst_degs <- asw_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder for plot
#asw_vst_degs <- asw_vst_degs[,c(1:22,24,27,28,32,37,38,41,42:79,23,25,26,29,30,31,33,34,35,36,39,40)]

##get tissue label info
sample_to_location <- data.table(data.frame(colData(asw_dds_para)[,c("Parasitism_status", "Weevil_Location", "sample_name")]))
sample_to_location <- sample_to_location %>% remove_rownames %>% column_to_rownames(var="sample_name")

location_colours <- list(Parasitism_status=c(parasitized="#FEAF77FF", undetected="#F1605DFF"), Weevil_Location = c(Dunedin="#3B0F70FF", Ruakura="#B63679FF"))
##plot
##not clustered by sample
pheatmap(asw_vst_degs, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE,
         annotation_col=sample_to_location, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))

