library(tximport)
library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
##factors and design
asw_dds$para <- factor(paste(asw_dds$Parasitism_status))
asw_dds$pc1_sign <- factor(paste(asw_dds$PC1_sign))
asw_dds_para <- copy(asw_dds)
##only Ru samples
asw_dds_para <- asw_dds_para[,asw_dds_para$Weevil_Location == "Dunedin"]
##control for PC1
design(asw_dds_para) <- ~pc1_sign+para
##run deseq2
asw_dds_para <- DESeq(asw_dds_para)
saveRDS(asw_dds_para, "output/deseq2/asw_dual/WT_parasitism_dunedin/parasitism_WT.rds")

asw_dds_para <- readRDS("output/deseq2/asw_dual/WT_parasitism_dunedin/parasitism_WT.rds")
##results
res_group <- results(asw_dds_para, lfcThreshold = 1, alpha = 0.1)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.1)
fwrite(ordered_sig_res_group_table, "output/deseq2/asw_dual/WT_parasitism_dunedin/sig_degs.csv")
fwrite(ordered_res_group_table, "output/deseq2/asw_dual/WT_parasitism_dunedin/res_group.csv")

##just ru and just dun find similar number DEGs as analysis controlling for location

##volcano plot
EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.1, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))


##merge with trinotate annots
trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw_dual/WT_parasitism_dunedin/sig_annots.csv")
