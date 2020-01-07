library("data.table")
library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("EnhancedVolcano")

dds_group <- readRDS("output/deseq2/dds_group.rds")

resultsNames(dds_group)

dun_res_group <- results(dds_group, contrast = c("group", "Dunedin_E", "Dunedin_N"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
dun_ordered_res_group <- dun_res_group[order(dun_res_group$padj),]
##Make data table and write to output
dun_ordered_res_group_table <- data.table(data.frame(dun_ordered_res_group), keep.rownames = TRUE)
dun_ordered_sig_res_group_table <- subset(dun_ordered_res_group_table, padj < 0.05)

##write full and sig res files
fwrite(dun_ordered_res_group_table, "output/deseq2/dunedin/res_group.csv")
fwrite(dun_ordered_sig_res_group_table, "output/deseq2/dunedin/dunedin_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)

##Sub in any gene of interest to plot counts  
plotCounts(dds_group, "TRINITY_DN2421_c0_g1", intgroup = c("group"), main="")
##volcano plot
EnhancedVolcano(dun_ordered_res_group_table, x="log2FoldChange", y="padj", lab="", pointSize = 3)

##read in annotated transcriptome
trinotate_report <- fread("data/asw_transcriptome/trinotate_annotation_report.txt")
setnames(dun_ordered_sig_res_group_table, old=c("rn"), new=c("#gene_id"))
##merge list of sig genes with annotations
dun_sig_w_annots <- merge (dun_ordered_sig_res_group_table, trinotate_report, by.x="#gene_id", by.y="#gene_id")
##save file - in excel edit duplicated gene ids (where one DEG had multiple annotations for each isoform)
fwrite(dun_sig_w_annots, "output/deseq2/dunedin/dun_sig_genes_with_annots.csv")
