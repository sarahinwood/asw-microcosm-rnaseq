library("tximport")
library("data.table")
library("DESeq2")

##~200 sig terms - why? Should be a comp between para ASW from 2 different locations..

dds_group_parasitism <- readRDS("output/old_para_transcriptome_deseq/parasitism_status/dds_group_parasitism.rds")
resultsNames(dds_group_parasitism)

para_res_group <- results(dds_group_parasitism, contrast = c("group", "Ruakura_Yes", "Dunedin_Yes"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
para_ordered_res_group <- para_res_group[order(para_res_group$padj),]
##Make data table and write to output
para_ordered_res_group_table <- data.table(data.frame(para_ordered_res_group), keep.rownames = TRUE)
para_ordered_sig_res_group_table <- subset(para_ordered_res_group_table, padj < 0.05)
para_sig_annots <- merge(para_ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE, all.y=FALSE)
##write full and sig res files
fwrite(para_ordered_res_group_table, "output/deseq2/parasitism_status/para_res_group.csv")
fwrite(para_sig_annots, "output/deseq2/parasitism_status/para_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)


plotCounts(dds_group_parasitism, "TRINITY_DN46809_c0_g1", intgroup = c("group"))



