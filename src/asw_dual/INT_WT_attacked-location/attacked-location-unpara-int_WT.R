library(tximport)
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
##relevel factors
asw_dds$attacked <- relevel(asw_dds$attacked, ref="Y")
asw_dds$location <- relevel(asw_dds$location, ref="Dunedin")
##remove para samples
asw_dds <- asw_dds[,asw_dds$para=="undetected"]


asw_dds_att_loc <- copy(asw_dds)
design(asw_dds_att_loc) <- ~pc1_sign+location+attacked+location:attacked
##run deseq2
asw_dds_att_loc <- DESeq(asw_dds_att_loc)
saveRDS(asw_dds_att_loc, "output/deseq2/asw_dual/INT_WT_attacked-location-unpara/parasitism-location-int_WT.rds")

asw_dds_att_loc <- readRDS("output/deseq2/asw_dual/INT_WT_attacked-location-unpara/parasitism-location-int_WT.rds")
##results
res_group <- results(asw_dds_att_loc, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/asw_dual/INT_WT_attacked-location-unpara/sig_degs.csv")
fwrite(ordered_res_group_table, "output/deseq2/asw_dual/INT_WT_attacked-location-unpara/res_group.csv")

asw_dds_att_loc$group <- factor(paste(asw_dds_att_loc$Attacked, asw_dds_att_loc$location, sep=" "))

plotCounts(asw_dds_att_loc, "ASW_TRINITY_DN18602_c0_g1", intgroup=("group"))

trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw_dual/INT_WT_attacked-location-unpara/sig_annots.csv")

#############
## heatmap ##
#############
all_sig_annots <- sig_annots
##vst transform
asw_vst <- varianceStabilizingTransformation(asw_dds_att_loc, blind=TRUE)
asw_vst_assay_dt <- data.table(assay(asw_vst), keep.rownames=TRUE)
##subset for DEGs
asw_vst_degs <- subset(asw_vst_assay_dt, rn %in% all_sig_annots$rn)
asw_vst_degs$rn <- tstrsplit(asw_vst_degs$rn, "ASW_", keep=c(2))
##turn first row back to row name
asw_vst_degs <- asw_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder for plot
asw_vst_degs <- asw_vst_degs[,c(1:22,24,27,28,32,37,38,41,23,25,26,29,30,31,33,34,35,36,39,40, 42:79)]

##get tissue label info
sample_to_location <- data.table(data.frame(colData(asw_dds_att_loc)[,c("location", "Attacked", "sample_name")]))
sample_to_location <- sample_to_location %>% remove_rownames %>% column_to_rownames(var="sample_name")

location_colours <- list(location = c(Dunedin="#3B0F70FF", Ruakura="#B63679FF"), Attacked=c(Y="#FEAF77FF", NO="#F1605DFF"))
##plot
##not clustered by sample
pheatmap(asw_vst_degs, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
         annotation_col=sample_to_location, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))

