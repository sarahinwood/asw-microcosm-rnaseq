library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(viridis)
library(pheatmap)
library(tidyverse)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
##factors and design
asw_dds$attacked <- factor(paste(asw_dds$Attacked))
asw_dds$para <- factor(paste(asw_dds$Parasitism_status))
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$pc1_sign <- factor(paste(asw_dds$PC1_sign))
##set reference level
asw_dds$attacked <- relevel(asw_dds$attacked, ref="Y")

asw_dds_attack <- copy(asw_dds)
design(asw_dds_attack) <- ~pc1_sign+para+location+attacked
##run deseq2
asw_dds_attack <- DESeq(asw_dds_attack)
saveRDS(asw_dds_attack, "output/deseq2/asw_dual/WT_attacked/attacked_WT.rds")

asw_dds_attack <- readRDS("output/deseq2/asw_dual/WT_attacked/attacked_WT.rds")
##results
res_group <- results(asw_dds_attack, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/asw_dual/WT_attacked/sig_degs.csv")
fwrite(ordered_res_group_table, "output/deseq2/asw_dual/WT_attacked/res_group.csv")

trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw_dual/WT_attacked/sig_annots.csv")

plotCounts(asw_dds_attack, "ASW_TRINITY_DN13344_c1_g1", intgroup="attacked")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.05, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

#############
## heatmap ##
#############
all_sig_annots <- fread("output/deseq2/asw_dual/WT_attacked/sig_blast_annots.csv")
##vst transform
asw_vst <- varianceStabilizingTransformation(asw_dds_attack, blind=TRUE)
asw_vst_assay_dt <- data.table(assay(asw_vst), keep.rownames=TRUE)
##subset for DEGs
asw_vst_degs <- subset(asw_vst_assay_dt, rn %in% all_sig_annots$rn)
asw_vst_degs$rn <- tstrsplit(asw_vst_degs$rn, "ASW_", keep=c(2))
##turn first row back to row name
asw_vst_degs <- asw_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder for plot
asw_vst_degs <- asw_vst_degs[,c(1:22,24,27,28,32,37,38,41,42:79,23,25,26,29,30,31,33,34,35,36,39,40)]

##get tissue label info
sample_to_location <- data.table(data.frame(colData(asw_dds_attack)[,c("Attacked", "sample_name")]))
sample_to_location <- sample_to_location %>% remove_rownames %>% column_to_rownames(var="sample_name")

location_colours <- list(Attacked=c(Y="#FEAF77FF", NO="#F1605DFF"))
##plot
##not clustered by sample
pheatmap(asw_vst_degs, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=TRUE,
         annotation_col=sample_to_location, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
