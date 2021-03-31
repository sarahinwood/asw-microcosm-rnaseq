library(tximport)
library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
##factors and design
asw_dds$attacked <- factor(paste(asw_dds$Attacked))
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$pc1_sign <- factor(paste(asw_dds$PC1_sign))
asw_dds_attack <- copy(asw_dds)
design(asw_dds_attack) <- ~location+pc1_sign+attacked
##run deseq2
asw_dds_attack <- DESeq(asw_dds_attack)
saveRDS(asw_dds_attack, "output/deseq2/asw_dual/WT_attacked/attacked_WT.rds")

asw_dds_attack <- readRDS("output/deseq2/asw_dual/WT_attacked/attacked_WT.rds")
##results
res_group <- results(asw_dds_attack, lfcThreshold = 1, alpha = 0.1)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.1)
fwrite(ordered_sig_res_group_table, "output/deseq2/asw_dual/attacked_WT/sig_degs.csv")
fwrite(ordered_res_group_table, "output/deseq2/asw_dual/WT_attacked/res_group.csv")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.05, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw_dual/WT_attacked/sig_annots.csv")
