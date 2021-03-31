library(data.table)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(EnhancedVolcano)

##Pairwise between cleaned and uncleaned
##only 13 DEGs - cleaning has not had huge result on samples

dds <- readRDS("output/deseq2/asw/asw_dds.rds")

##create dds object for cleaned analysis
dds_cleaned <- copy(dds)
dds_cleaned$group <- factor(paste(dds$Cleaned))

##add group to design
design(dds_cleaned) <- ~group
##run deseq2 and generate results
dds_cleaned <- DESeq(dds_cleaned)
##save dds_group
saveRDS(dds_cleaned, file = "output/deseq2/asw/WT_cleaned/asw_dds_cleaned.rds")

resultsNames(dds_cleaned)

##Make table of results for exposed vs control heads
res_group <- results(dds_cleaned, contrast = c("group", "Yes", "No"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
##write tables
fwrite(ordered_res_group_table, "output/deseq2/asw/WT_cleaned/full_res.csv")
fwrite(ordered_sig_res_group_table, "output/deseq2/asw/WT_cleaned/sig_degs.csv", col.names = TRUE, row.names = FALSE)

###merge with annots
longest_trinotate_report <- fread("data/asw_edited_transcript_ids/trinotate_longest_isoform.csv")
sig_degs_annots <- merge(ordered_sig_res_group_table, longest_trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_degs_annots, "output/deseq2/asw/WT_cleaned/sig_degs_annots.csv")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", pointSize = 3, pCutoff=0.05)

dds_cleaned$loc_para <- factor(paste(dds$Weevil_Location, dds$Parasitism_status, sep="_"))
##can add in cleaned to check DE isn't a result of cleaned
plotCounts(dds_cleaned, "ASW_TRINITY_DN29024_c0_g1", intgroup = c("group", "loc_para"))
