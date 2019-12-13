library(data.table)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(EnhancedVolcano)
library(VennDiagram)

dds_group <- readRDS("output/concat_deseq2/dds_group.rds")

resultsNames(dds_group)

dun_res_group <- results(dds_group, contrast = c("group", "Dunedin_E", "Dunedin_N"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
dun_ordered_res_group <- dun_res_group[order(dun_res_group$padj),]
##Make data table and write to output
dun_ordered_res_group_table <- data.table(data.frame(dun_ordered_res_group), keep.rownames = TRUE)
dun_ordered_sig_res_group_table <- subset(dun_ordered_res_group_table, padj < 0.05)