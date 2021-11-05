library(data.table)
library(tidyverse)
library(DESeq2)

mh_dds <- readRDS("output/deseq2/mh_dual/mh_dual_dds.rds")
trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings=".")

##factors and design
mh_dds$pc1_sign <- factor(paste(mh_dds$PC1_sign))
mh_dds$para <- factor(mh_dds$Parasitism_status)
mh_dds$viral <- factor(mh_dds$Viral_expressed)

mh_dds_viral <- copy(mh_dds)
##remove para samples
mh_dds_viral <- mh_dds_viral[,mh_dds_viral$para=="undetected"]
design(mh_dds_viral) <- ~pc1_sign+viral
##run deseq2
mh_dds_viral <- DESeq(mh_dds_viral)

###########
##results##
###########
resultsNames(mh_dds_viral)
viral_res_group <- results(mh_dds_viral, contrast=c("viral", "Yes", "No"), lfcThreshold = 1, alpha = 0.05)
##3 DEGs
summary(viral_res_group)
##Order based of padj
viral_ordered_res_group <- viral_res_group[order(viral_res_group$padj),]
##Make data table and write to output
viral_ordered_res_group_table <- data.table(data.frame(viral_ordered_res_group), keep.rownames = TRUE)
viral_ordered_sig_res_group_table <- subset(viral_ordered_res_group_table, padj < 0.05)
##DEG annotations
viral_sig_annots <- merge(viral_ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
##no annotations
fwrite(viral_sig_annots, "output/deseq2/mh_dual/unpara_viral_WT/sig_annots.csv")
plotCounts(mh_dds_viral, "MH_TRINITY_DN85248_c0_g1", intgroup="viral")
