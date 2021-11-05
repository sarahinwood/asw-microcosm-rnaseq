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
asw_dds$pc1_sign <- factor(paste(asw_dds$PC1_sign))
##relevel parasitism
asw_dds$para <- relevel(asw_dds$para, ref="undetected")
##only Invermay
asw_dds_inv <- asw_dds[,asw_dds$Weevil_Location=="Dunedin"]

asw_dds_para <- copy(asw_dds_inv)
design(asw_dds_para) <- ~pc1_sign+para
##run deseq2
asw_dds_para <- DESeq(asw_dds_para)

##results
res_group <- results(asw_dds_para, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)

##overlap with all para asw analysis
all_sig <- fread("output/deseq2/asw_dual/WT_parasitism/sig_blast_annots.csv")
inv_only <- setdiff(ordered_sig_res_group_table$rn, all_sig$rn)
all_only <- setdiff(all_sig$rn, ordered_sig_res_group_table$rn)
all_only_annots <- subset(all_sig, rn %in% all_only)

plotCounts(asw_dds_para, "ASW_TRINITY_DN6332_c1_g1", intgroup="para")
##volcano plot
EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.1, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

##merge with trinotate annots
trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
sig_annots_inv_only <- subset(sig_annots, rn %in% inv_only)
