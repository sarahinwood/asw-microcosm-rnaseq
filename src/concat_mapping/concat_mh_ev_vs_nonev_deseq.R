library(data.table)
library(DESeq2)

dds_concat_group_mh <- readRDS("output/concat_deseq2/dds_concat_behaviour_mh.rds")
resultsNames(dds_concat_group_mh)

res_group <- results(dds_concat_group_mh, contrast = c("group", "E", "N"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
