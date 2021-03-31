library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(VennDiagram)

##Pairwise between unpara-viral-expressing and unpara-no-viral

dds <- readRDS("output/deseq2/asw/asw_dds.rds")

##create dds object for parasitism analysis
dds_para_viral <- copy(dds)
dds_para_viral$group <- factor(paste(dds$Parasitism_status, dds$Viral_expressed, sep="_"))

##add group to design
design(dds_para_viral) <- ~group
##run deseq2 and generate results
dds_para_viral <- DESeq(dds_para_viral)
##save dds_group
saveRDS(dds_para_viral, file = "output/deseq2/asw/WT_unpara_viral/asw_dds_para_viral.rds")

resultsNames(dds_para_viral)

##Make table of results for exposed vs control heads
res_group <- results(dds_para_viral, contrast = c("group", "undetected_Yes", "undetected_No"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
##write tables
fwrite(ordered_res_group_table, "output/deseq2/asw/WT_unpara_viral/full_res.csv")
fwrite(ordered_sig_res_group_table, "output/deseq2/asw/WT_unpara_viral/sig_degs.csv", col.names = TRUE, row.names = FALSE)

###merge with annots
longest_trinotate_report <- fread("data/asw_edited_transcript_ids/trinotate_longest_isoform.csv")
sig_degs_annots <- merge(ordered_sig_res_group_table, longest_trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_degs_annots, "output/deseq2/asw/WT_unpara_viral/sig_degs_annots.csv")

EnhancedVolcano(ordered_res_group_table, x="log2FoldChange", y="padj", lab="", pointSize = 3, pCutoff=0.05)

##can add in parasitism to check DE isn't a result of parasitism
plotCounts(dds_para_viral, "asw_TRINITY_DN21057_c0_g1", intgroup = c("group"))

para_sig <- fread("output/deseq2/asw/WT_parasitism/sig_degs.csv")
intersect(ordered_sig_res_group_table$rn, para_sig$rn)

vd <- venn.diagram(x = list("parasitism DEGs"=para_sig$rn, "unpara viral DEGs"=ordered_sig_res_group_table$rn), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1,)
grid.newpage()
grid.draw(vd)
