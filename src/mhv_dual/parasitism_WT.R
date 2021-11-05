library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(dplyr)
library(pheatmap)

mhv_dds <- readRDS("output/deseq2/mhv_dual/mhv_dual_dds.rds")
##factors and design
mhv_dds$para <- factor(paste(mhv_dds$Parasitism_status))
mhv_dds$para <- relevel(mhv_dds$para, ref="undetected")
mhv_dds$location <- factor(paste(mhv_dds$Weevil_Location))
mhv_dds$pc1_sign <- factor(paste(mhv_dds$PC1_sign))
mhv_dds_para <- copy(mhv_dds)
design(mhv_dds_para) <- ~location+pc1_sign+para
##run deseq2
mhv_dds_para <- DESeq(mhv_dds_para)
saveRDS(mhv_dds_para, "output/deseq2/mhv_dual/parasitism_WT/parasitism_WT.rds")

mhv_dds_para <- readRDS("output/deseq2/mhv_dual/parasitism_WT/parasitism_WT.rds")
##results
res_group <- results(mhv_dds_para, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/mhv_dual/parasitism_WT/sig_degs.csv")
fwrite(ordered_res_group_table, "output/deseq2/mhv_dual/parasitism_WT/res_group.csv")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", title="",
                subtitle="", pointSize = 1.5, pCutoff = 0.05, colAlpha=0.2,
                col=c("#FDE725FF", "#21908CFF", "grey20", "#440154FF"))

trinotate_report <- fread("data/mhv-mh-combined-transcriptome/output/mh_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
##all sig with trinotate annots
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/mhv_dual/parasitism_WT/sig_annots.csv")
##merge with MhV annots
mhv_genes <- fread("data/mh-rnaseq/output/blast/viral_genes/viral_genes_best_hits.csv")
mhv_genes$full_id <- paste("MH", mhv_genes$Trinity_ID, sep="_")
mhv_genes$full_id <- tstrsplit(mhv_genes$full_id, "_i", keep=c(1))
sig_mhv_annots <- merge(sig_annots, mhv_genes, by.x="rn", by.y="full_id")
fwrite(sig_mhv_annots, "output/deseq2/mhv_dual/parasitism_WT/sig_MhV_annots.csv")

sig_annots$edited_rn <- tstrsplit(sig_annots$rn, "MH_", keep=c(2))
mhv_genes_de <- subset(mhv_genes, mhv_genes$`#gene_id` %in% sig_annots$edited_rn)

#############
## heatmap ##
#############

##vst transform
mhv_vst <- varianceStabilizingTransformation(mhv_dds_para, blind=TRUE)
mhv_vst_assay_dt <- data.table(assay(mhv_vst), keep.rownames=TRUE)
##subset for DEGs
mhv_vst_degs <- subset(mhv_vst_assay_dt, rn %in% sig_annots$rn)
mhv_vst_degs$rn <- tstrsplit(mhv_vst_degs$rn, "mhv_", keep=c(2))
##turn first row back to row name
mhv_vst_degs <- mhv_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")

##get tissue label info
sample_to_location <- data.table(data.frame(colData(mhv_dds_para)[,c("Parasitism_status", "location", "sample_name")]))
sample_to_location <- sample_to_location %>% remove_rownames %>% column_to_rownames(var="sample_name")

location_colours <- list(Parasitism_status = c(parasitized="#FEAF77FF", undetected="#F1605DFF"), location = c(Dunedin="#3B0F70FF", Ruakura="#B63679FF"))
##plot
##clustered by sample
pheatmap(mhv_vst_degs, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE,
         annotation_col=sample_to_location, annotation_colors=location_colours, annotation_names_col=FALSE,
         show_colnames = TRUE, border_color=NA, color=viridis(50))

##plot counts against PCR target for five samples with viral expression
sample_data <- fread("data/sample_table.csv")
unpara <- subset(sample_data, sample_data$Parasitism_status=="undetected")
unpara_viral <- subset(unpara, unpara$Viral_expressed=="Yes")

mh_dds_para <- readRDS("output/deseq2/mh_dual/parasitism_WT/parasitism_WT.rds")
PCR_counts <- data.table(plotCounts(mh_dds_para, "MH_TRINITY_DN481_c1_g1", intgroup="para", returnData=TRUE), keep.rownames = TRUE)
unpara_viralcounts <- subset(PCR_counts, rn %in% unpara_viral$sample_name)

