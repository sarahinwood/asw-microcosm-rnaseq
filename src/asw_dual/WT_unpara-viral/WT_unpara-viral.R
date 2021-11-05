library(data.table)
library(tidyverse)
library(DESeq2)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings=".")

##factors and design
asw_dds$pc1_sign <- factor(paste(asw_dds$PC1_sign))
asw_dds$para <- factor(asw_dds$Parasitism_status)
asw_dds$viral <- factor(asw_dds$Viral_expressed)

asw_dds_viral <- copy(asw_dds)
##remove para samples
asw_dds_viral <- asw_dds_viral[,asw_dds_viral$para=="undetected"]
design(asw_dds_viral) <- ~pc1_sign+viral
##run deseq2
asw_dds_viral <- DESeq(asw_dds_viral)

###########
##results##
###########
resultsNames(asw_dds_viral)
viral_res_group <- results(asw_dds_viral, contrast=c("viral", "Yes", "No"), lfcThreshold = 1, alpha = 0.05)
##3 DEGs
summary(viral_res_group)
##Order based of padj
viral_ordered_res_group <- viral_res_group[order(viral_res_group$padj),]
##Make data table and write to output
viral_ordered_res_group_table <- data.table(data.frame(viral_ordered_res_group), keep.rownames = TRUE)
viral_ordered_sig_res_group_table <- subset(viral_ordered_res_group_table, padj < 0.05)
##DEG annotations
viral_sig_annots <- merge(viral_ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
##no annotations

plotCounts(asw_dds_viral, "ASW_TRINITY_DN6075_c0_g1", intgroup="viral")
