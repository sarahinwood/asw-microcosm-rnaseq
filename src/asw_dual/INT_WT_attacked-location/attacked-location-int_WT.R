library(tximport)
library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
##factors and design
asw_dds$attacked <- factor(paste(asw_dds$Attacked))
asw_dds$para <- factor(paste(asw_dds$Parasitism_status))
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$pc1_sign <- factor(paste(asw_dds$PC1_sign))
##relevel factors
asw_dds$attacked <- relevel(asw_dds$attacked, ref="Y")
asw_dds$location <- relevel(asw_dds$location, ref="Dunedin")

asw_dds_att_loc <- copy(asw_dds)
design(asw_dds_att_loc) <- ~pc1_sign+para+location+attacked+location:attacked
##run deseq2
asw_dds_att_loc <- DESeq(asw_dds_att_loc)
saveRDS(asw_dds_att_loc, "output/deseq2/asw_dual/INT_WT_attacked-location/parasitism-location-int_WT.rds")

asw_dds_att_loc <- readRDS("output/deseq2/asw_dual/INT_WT_attacked-location/parasitism-location-int_WT.rds")
##results
res_group <- results(asw_dds_att_loc, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/asw_dual/INT_WT_attacked-location/sig_degs.csv")
fwrite(ordered_res_group_table, "output/deseq2/asw_dual/INT_WT_attacked-location/res_group.csv")

asw_dds_att_loc$group <- factor(paste(asw_dds_att_loc$Attack_status, asw_dds_att_loc$location, sep=" "))
plotCounts(asw_dds_att_loc, "ASW_TRINITY_DN17902_c1_g2", intgroup=("group"))

trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw_dual/INT_WT_attacked-location/sig_annots.csv")
