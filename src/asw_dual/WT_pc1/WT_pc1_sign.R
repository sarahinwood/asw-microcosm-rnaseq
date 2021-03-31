library(tximport)
library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")

##factors and design
asw_dds$para <- factor(paste(asw_dds$Parasitism_status))
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$pc1_sign <- factor(paste(asw_dds$PC1_sign))
asw_dds_pc1 <- copy(asw_dds)
design(asw_dds_pc1) <- ~para+location+pc1_sign
##run deseq2
asw_dds_pc1 <- DESeq(asw_dds_pc1)

##results
res_group <- results(asw_dds_pc1, contrast = c("pc1_sign", "negative", "positive"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
fwrite(ordered_res_group_table, "output/deseq2/asw_dual/WT_pc1/res_group.csv")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.05, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/asw_dual/WT_pc1/sig_degs.csv")

trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw_dual/WT_pc1/sig_annots.csv")
