library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(dplyr)
library(viridis)
library(pheatmap)
library(tidyverse)

mh_dds <- readRDS("output/deseq2/mh_dual/mh_dual_dds.rds")
##factors and design
mh_dds$para <- factor(paste(mh_dds$Parasitism_status))
mh_dds$para <- relevel(mh_dds$para, ref="undetected")
mh_dds$location <- factor(paste(mh_dds$Weevil_Location))
mh_dds$pc1_sign <- factor(paste(mh_dds$PC1_sign))
mh_dds_para <- copy(mh_dds)
design(mh_dds_para) <- ~location+pc1_sign+para
##run deseq2
mh_dds_para <- DESeq(mh_dds_para)
saveRDS(mh_dds_para, "output/deseq2/mh_dual/parasitism_WT/parasitism_WT.rds")

mh_dds_para <- readRDS("output/deseq2/mh_dual/parasitism_WT/parasitism_WT.rds")
##results
res_group <- results(mh_dds_para, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/mh_dual/parasitism_WT/sig_degs.csv")
fwrite(ordered_res_group_table, "output/deseq2/mh_dual/parasitism_WT/res_group.csv")

plotCounts(mh_dds_para, "MH_TRINITY_DN1562_c1_g1", intgroup="para")
EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.05, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/mh_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
##all sig with trinotate annots
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/mh_dual/parasitism_WT/sig_annots.csv")

venom_degs <- fread("/Volumes/archive/deardenlab/sarahinwood/mh_projects/mh-rnaseq/output/deseq2/tissue_itWT_LRT/venom/venom_sp_LRT_annots.csv")
up_venom <- subset(venom_degs, log2FoldChange<0)
up_venom_sp <- subset(up_venom, SignalP=="Y")

##are crawford venom genes expressed by unpara_viral ASW?
crawford_degs <- fread("/Volumes/archive/deardenlab/sarahinwood/mh_projects/mh-rnaseq/output/blast/crawford_venom/crawford_best_hits.csv")
##ASW. VG2, VG3 and VG4 are not upregulated in parasitized ASW
mh_dds_para$group <- factor(paste(mh_dds_para$Viral_expressed, mh_dds_para$Parasitism_status, sep="_"))
plotCounts(mh_dds_para, "MH_TRINITY_DN200_c0_g1", intgroup="group", main = "VG1")
plotCounts(mh_dds_para, "MH_TRINITY_DN2387_c0_g1", intgroup="group", main = "VG5")
plotCounts(mh_dds_para, "MH_TRINITY_DN1562_c1_g1", intgroup="group", main = "VG6")
plotCounts(mh_dds_para, "MH_TRINITY_DN5831_c0_g1", intgroup="group", main = "VG8")
plotCounts(mh_dds_para, "MH_TRINITY_DN5069_c0_g1", intgroup="group", main = "VG10")

sig_annots$edited_rn <- tstrsplit(sig_annots$rn, "MH_", keep=c(2))
venom_sp_degs <- intersect(up_venom_sp$rn, sig_annots$edited_rn)
##includes 47 of the venom DEGs with SPs
venom_sp_degs_annots <- subset(venom_degs, rn %in% venom_sp_degs)

#############
## heatmap ##
#############
##vst transform
mh_vst <- varianceStabilizingTransformation(mh_dds_para, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
##subset for DEGs
mh_vst_assay_dt$rn <- tstrsplit(mh_vst_assay_dt$rn, "MH_", keep=c(2))
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% venom_sp_degs)
##turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
##reorder for plot
#mh_vst_degs <- mh_vst_degs[,c(1:22,24,27,28,32,37,38,41,42:79,23,25,26,29,30,31,33,34,35,36,39,40)]

##get tissue label info
sample_to_location <- data.table(data.frame(colData(mh_dds_para)[,c("Parasitism_status", "sample_name")]))
sample_to_location <- sample_to_location %>% remove_rownames %>% column_to_rownames(var="sample_name")

location_colours <- list(Parasitism_status=c(parasitized="#FEAF77FF", undetected="#F1605DFF"))
##plot
##not clustered by sample
pheatmap(mh_vst_degs, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE,
         annotation_col=sample_to_location, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
