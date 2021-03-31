library(tximport)
library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
##factors and design
asw_dds$attack <- factor(paste(asw_dds$Attack_status))
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$pc1_sign <- factor(paste(asw_dds$PC1_sign))
asw_dds_att_stat_loc <- copy(asw_dds)

##Look for genes DE between locations, controlling for attack status
design(asw_dds_att_stat_loc) <- ~pc1_sign+attack+location+attack:location
asw_dds_att_stat_loc <- DESeq(asw_dds_att_stat_loc, test='LRT', reduced=~pc1_sign+attack+location+attack)
saveRDS(asw_dds_att_stat_loc, "output/deseq2/asw_dual/INT_LRT_attack-status-location/att-stat-location-int_WT.rds")

asw_dds_att_stat_loc <- readRDS("output/deseq2/asw_dual/INT_LRT_attack-status-location/att-stat-location-int_WT.rds")
##results
res_group <- results(asw_dds_att_stat_loc, alpha = 0.1)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.1)
fwrite(ordered_sig_res_group_table, "output/deseq2/asw_dual/INT_LRT_attack-status-location/sig_degs.csv")
fwrite(ordered_res_group_table, "output/deseq2/asw_dual/INT_LRT_attack-status-location/res_group.csv")

trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw_dual/INT_LRT_attack-status-location/sig_annots.csv")
