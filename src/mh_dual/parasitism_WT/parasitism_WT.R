library(tximport)
library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(dplyr)

mh_dds <- readRDS("output/deseq2/mh_dual/mh_dual_dds.rds")
##factors and design
mh_dds$para <- factor(paste(mh_dds$Parasitism_status))
mh_dds$location <- factor(paste(mh_dds$Weevil_Location))
mh_dds$pc1_sign <- factor(paste(mh_dds$PC1_sign))
mh_dds_para <- copy(mh_dds)
design(mh_dds_para) <- ~location+pc1_sign+para
##run deseq2
mh_dds_para <- DESeq(mh_dds_para)
saveRDS(mh_dds_para, "output/deseq2/mh_dual/parasitism_WT/parasitism_WT.rds")

mh_dds_para <- readRDS("output/deseq2/mh_dual/parasitism_WT/parasitism_WT.rds")
##results
res_group <- results(mh_dds_para, contrast=c("para", "parasitised", "undetected"), lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/mh_dual/parasitism_WT/sig_degs.csv")
fwrite(ordered_res_group_table, "output/deseq2/mh_dual/parasitism_WT/res_group.csv")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.05, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/mh_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
viral_recip <- fread("/Volumes/archive/deardenlab/sarahinwood/mh_projects/mh-transcriptome/output/recip_blast/nr_blastx/viral_transcripts_annots.csv")
viral_recip$edited_ids <- paste("MH", viral_recip$transcript_id, sep = "_")
viral_recip$edited_ids <- tstrsplit(viral_recip$edited_ids, "_i", keep=c(1))
##all sig with trinotate annots
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/mh_dual/parasitism_WT/sig_annots.csv")
##sig degs with viral blastx
sig_annots_recip <- merge(sig_annots, viral_recip, by.x="rn", by.y="edited_ids", all.x=TRUE)
fwrite(sig_annots_recip, "output/deseq2/mh_dual/parasitism_WT/sig_recip_viral_annots.csv")

virus_x <- data.table(dplyr::filter(sig_annots_recip, grepl('Viruses', sprot_Top_BLASTX_hit)))
virus_y <- data.table(dplyr::filter(sig_annots_recip, grepl('virus', annotation)))
viral_degs <- full_join(virus_x, virus_y)
fwrite(viral_degs, "output/deseq2/mh_dual/parasitism_WT/sig_viral_degs.csv")

