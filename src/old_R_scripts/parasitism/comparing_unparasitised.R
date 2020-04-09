library("tximport")
library("data.table")
library("DESeq2")

dds_group_parasitism <- readRDS("output/old_para_transcriptome_deseq/parasitism_status/dds_group_parasitism.rds")
resultsNames(dds_group_parasitism)

unpara_res_group <- results(dds_group_parasitism, contrast = c("group", "Ruakura_No", "Dunedin_No"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
unpara_ordered_res_group <- unpara_res_group[order(unpara_res_group$padj),]
##Make data table and write to output
unpara_ordered_res_group_table <- data.table(data.frame(unpara_ordered_res_group), keep.rownames = TRUE)
unpara_ordered_sig_res_group_table <- subset(unpara_ordered_res_group_table, padj < 0.05)
unpara_sig_annots <- merge(unpara_ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE, all.y=FALSE)
##write full and sig res files
fwrite(unpara_ordered_res_group_table, "output/deseq2/parasitism_status/unpara_res_group.csv")
fwrite(unpara_sig_annots, "output/deseq2/parasitism_status/unpara_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)


plotCounts(dds_group_parasitism, "TRINITY_DN18865_c0_g1", intgroup = c("group"))



