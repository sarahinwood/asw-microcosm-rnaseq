library(data.table)
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(viridis)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
##factors and design
asw_dds$para <- factor(paste(asw_dds$Parasitism_status))
asw_dds$pc1_sign <- factor(paste(asw_dds$PC1_sign))
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
##set reference level
asw_dds$location <- relevel(asw_dds$location, ref="Dunedin")

asw_dds_location <- copy(asw_dds)
design(asw_dds_location) <- ~pc1_sign+para+location
##run deseq2
asw_dds_location <- DESeq(asw_dds_location)
saveRDS(asw_dds_location, "output/deseq2/asw_dual/WT_location/location_WT.rds")

asw_dds_location <- readRDS("output/deseq2/asw_dual/WT_location/location_WT.rds")
##results
res_group <- results(asw_dds_location, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/asw_dual/WT_location/sig_degs.csv")
fwrite(ordered_res_group_table, "output/deseq2/asw_dual/WT_location/res_group.csv")

##DEG annotations
trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings=".")
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw_dual/WT_location/sig_annots.csv")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.1, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

plotCounts(asw_dds_location, "ASW_TRINITY_DN8550_c0_g1", intgroup="location")

##overlap with exposed-RNA-seq location analysis
exposed_location_degs <- fread("/Volumes/archive/deardenlab/sarahinwood/asw_projects/asw-exposed-rnaseq/output/deseq2/asw/all_location/sig_w_annots.csv")
ordered_sig_res_group_table$rn_edited <- tstrsplit(ordered_sig_res_group_table$rn, "ASW_", keep=c(2))
shared_degs <- intersect(exposed_location_degs$rn, ordered_sig_res_group_table$rn_edited)
exposed_sp <- setdiff(exposed_location_degs$rn, ordered_sig_res_group_table$rn_edited)
evasion_sp <- setdiff(ordered_sig_res_group_table$rn_edited, exposed_location_degs$rn)
##which degs in both?
sig_annots$edited_rn <- tstrsplit(sig_annots$rn, "ASW_", keep=c(2))
shared_sig_annots <- subset(sig_annots, edited_rn %in% shared_degs)

#############
## heatmap ##
#############
all_sig_annots <- fread("output/deseq2/asw_dual/WT_location/sig_blast_annots.csv")
##vst transform
asw_vst <- varianceStabilizingTransformation(asw_dds_location, blind=TRUE)
asw_vst_assay_dt <- data.table(assay(asw_vst), keep.rownames=TRUE)
##subset for DEGs
asw_vst_degs <- subset(asw_vst_assay_dt, rn %in% all_sig_annots$rn)
asw_vst_degs$rn <- tstrsplit(asw_vst_degs$rn, "ASW_", keep=c(2))
##turn first row back to row name
asw_vst_degs <- asw_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")

##get tissue label info
sample_to_location <- data.table(data.frame(colData(asw_dds_location)[,c("location", "sample_name")]))
sample_to_location <- sample_to_location %>% remove_rownames %>% column_to_rownames(var="sample_name")

location_colours <- list(location = c(Dunedin="#3B0F70FF", Ruakura="#B63679FF"))
##plot
##not clustered by sample
pheatmap(asw_vst_degs, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
         annotation_col=sample_to_location, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))

##also plot only annotated
only_annot <- subset(all_sig_annots, sprot_Top_BLASTX_hit!="")
asw_vst_degs_annot <- subset(asw_vst_assay_dt, rn %in% only_annot$rn)
asw_vst_degs_annot$rn <- tstrsplit(asw_vst_degs_annot$rn, "ASW_", keep=c(2))
##turn first row back to row name
asw_vst_degs_annot <- asw_vst_degs_annot %>% remove_rownames %>% column_to_rownames(var="rn")
##plot
pheatmap(asw_vst_degs_annot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
         annotation_col=sample_to_location, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
