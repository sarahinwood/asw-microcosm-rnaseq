library(data.table)
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(viridis)
library(VennDiagram)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings=".")

##factors and design
asw_dds$pc1_sign <- factor(paste(asw_dds$PC1_sign))
asw_dds$para <- factor(asw_dds$Parasitism_status)
asw_dds$group <- factor(paste(asw_dds$Weevil_Location, asw_dds$Attacked, sep="_"))

asw_dds_group <- copy(asw_dds)
design(asw_dds_group) <- ~pc1_sign+group
##run deseq2
asw_dds_group <- DESeq(asw_dds_group)
#saveRDS(asw_dds_group, "output/deseq2/asw_dual/WT_location-attack/location-attack_WT.rds")

###########
##results##
###########
#asw_dds_group <- readRDS("output/deseq2/asw_dual/WT_location-attack/location-attack_WT.rds")
resultsNames(asw_dds_group)
##Invermay
inv_res_group <- results(asw_dds_group, contrast=c("group", "Dunedin_NO", "Dunedin_Y"), lfcThreshold = 1, alpha = 0.05)
summary(inv_res_group)
##Order based of padj
inv_ordered_res_group <- inv_res_group[order(inv_res_group$padj),]
##Make data table and write to output
inv_ordered_res_group_table <- data.table(data.frame(inv_ordered_res_group), keep.rownames = TRUE)
inv_ordered_sig_res_group_table <- subset(inv_ordered_res_group_table, padj < 0.05)
fwrite(inv_ordered_sig_res_group_table, "output/deseq2/asw_dual/WT_location-attack/Invermay_sig_degs.csv")
fwrite(inv_ordered_res_group_table, "output/deseq2/asw_dual/WT_location-attack/Invermay_res_group.csv")
##DEG annotations
inv_sig_annots <- merge(inv_ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(inv_sig_annots, "output/deseq2/asw_dual/WT_location-attack/Invermay_sig_annots.csv")

##Ruakura
ru_res_group <- results(asw_dds_group, contrast=c("group", "Ruakura_NO", "Ruakura_Y"), lfcThreshold = 1, alpha = 0.05)
summary(ru_res_group)
##Order based of padj
ru_ordered_res_group <- ru_res_group[order(ru_res_group$padj),]
##Make data table and write to output
ru_ordered_res_group_table <- data.table(data.frame(ru_ordered_res_group), keep.rownames = TRUE)
ru_ordered_sig_res_group_table <- subset(ru_ordered_res_group_table, padj < 0.05)
fwrite(ru_ordered_sig_res_group_table, "output/deseq2/asw_dual/WT_location-attack/Ruakura_sig_degs.csv")
fwrite(ru_ordered_res_group_table, "output/deseq2/asw_dual/WT_location-attack/Ruakura_res_group.csv")
##DEG annotations
ru_sig_annots <- merge(ru_ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(ru_sig_annots, "output/deseq2/asw_dual/WT_location-attack/Ruakura_sig_annots.csv")

shared <- intersect(ru_sig_annots$rn, inv_sig_annots$rn)

##overlap with location:attack
att_loc <- fread("output/deseq2/asw_dual/INT_WT_attacked-location/sig_annots.csv")
vd <- venn.diagram(x = list("Interaction analysis"=att_loc$rn, "Ruakura Attack"=ru_sig_annots$rn, "Invermay Attack"=inv_sig_annots$rn), filename=NULL,
                   fill=c("#440154FF", "#21908CFF", "#FDE725FF"), alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd)

##plots
EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.1, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

plotCounts(asw_dds_group, "ASW_TRINITY_DN16678_c0_g1", intgroup="group")

#############
## heatmap ##
#############
##get tissue label info
sample_to_location <- data.table(data.frame(colData(asw_dds)[,c("Weevil_Location", "Attacked", "sample_name", "para")]))
sample_to_location <- sample_to_location %>% remove_rownames %>% column_to_rownames(var="sample_name")
location_colours <- list(Weevil_Location = c(Dunedin="#3B0F70FF", Ruakura="#B63679FF"), Attacked=c(Y="#FEAF77FF", NO="#F1605DFF"))

##vst transform
asw_vst <- varianceStabilizingTransformation(asw_dds_group, blind=TRUE)
asw_vst_assay_dt <- data.table(assay(asw_vst), keep.rownames=TRUE)

### Invermay
inv_sig_annots <- fread("output/deseq2/asw_dual/WT_location-attack/Invermay_sig_annots.csv")
##subset for DEGs
inv_asw_vst_degs <- subset(asw_vst_assay_dt, rn %in% inv_sig_annots$rn)
inv_asw_vst_degs$rn <- tstrsplit(inv_asw_vst_degs$rn, "ASW_", keep=c(2))
##turn first row back to row name
inv_asw_vst_degs <- inv_asw_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder for plot
inv_asw_vst_degs <- inv_asw_vst_degs[,c(1:22,24,27,28,32,37,38,41,23,25,26,29,30,31,33,34,35,36,39,40)]
##plot
##not clustered by sample
pheatmap(inv_asw_vst_degs, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=TRUE,
         annotation_col=sample_to_location, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))

### Ruakura
ru_sig_annots <- fread("output/deseq2/asw_dual/WT_location-attack/Ruakura_sig_annots.csv")
##subset for DEGs
ru_asw_vst_degs <- subset(asw_vst_assay_dt, rn %in% ru_sig_annots$rn)
ru_asw_vst_degs$rn <- tstrsplit(ru_asw_vst_degs$rn, "ASW_", keep=c(2))
##turn first row back to row name
ru_asw_vst_degs <- ru_asw_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder for plot
ru_asw_vst_degs <- ru_asw_vst_degs[,c(42:79)]
##plot
##not clustered by sample
pheatmap(ru_asw_vst_degs, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=TRUE,
         annotation_col=sample_to_location, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))

