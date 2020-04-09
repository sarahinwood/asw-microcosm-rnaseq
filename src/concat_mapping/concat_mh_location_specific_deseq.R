library(data.table)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(EnhancedVolcano)
library(VennDiagram)

dds_concat_group_mh <- readRDS("output/concat_deseq2/dds_concat_group_mh.rds")

resultsNames(dds_concat_group_mh)

###########
##RUAKURA##
###########

rua_res_group <- results(dds_concat_group_mh, contrast = c("group", "Ruakura_E", "Ruakura_N"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
rua_ordered_res_group <- rua_res_group[order(rua_res_group$padj),]
##Make data table and write to output
rua_ordered_res_group_table <- data.table(data.frame(rua_ordered_res_group), keep.rownames = TRUE)
rua_ordered_sig_res_group_table <- subset(rua_ordered_res_group_table, padj < 0.05)
rua_ordered_sig_res_group_table$rn <- tstrsplit(rua_ordered_sig_res_group_table$rn, "MH_", keep=c(2))

##write full and sig res files
fwrite(rua_ordered_res_group_table, "output/concat_deseq2/ruakura/res_group.csv")
fwrite(rua_ordered_sig_res_group_table, "output/concat_deseq2/ruakura/ruakura_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)

###########
##DUNEDIN##
###########

dun_res_group <- results(dds_concat_group_mh, contrast = c("group", "Dunedin_E", "Dunedin_N"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
dun_ordered_res_group <- dun_res_group[order(dun_res_group$padj),]
##Make data table and write to output
dun_ordered_res_group_table <- data.table(data.frame(dun_ordered_res_group), keep.rownames = TRUE)
dun_ordered_sig_res_group_table <- subset(dun_ordered_res_group_table, padj < 0.05)
dun_ordered_sig_res_group_table$rn <- tstrsplit(dun_ordered_sig_res_group_table$rn, "MH_", keep=c(2))

##write full and sig res files
fwrite(dun_ordered_res_group_table, "output/concat_deseq2/dunedin/res_group.csv")
fwrite(dun_ordered_sig_res_group_table, "output/concat_deseq2/dunedin/dunedin_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)

###############

plotCounts(dds_concat_group_mh, "MH_TRINITY_DN126_c0_g2", intgroup = c("group"), main = "")

##read in annotated transcriptome
trinotate_most_sig_hits <- fread("data/asw_mh_concat_transcriptome/mh_trinotate_best_hit_per_gene.csv")
setnames(rua_ordered_sig_res_group_table, old=c("rn"), new=c("#gene_id"))
##merge list of sig genes with annotations
rua_sig_w_annots <- merge (rua_ordered_sig_res_group_table, trinotate_most_sig_hits, by="#gene_id", all.x=TRUE)
##save file - no need to dedup anymore
fwrite(rua_sig_w_annots, "output/concat_deseq2/ruakura/rua_sig_genes_with_annots.csv")

not_concat_ru_degs <- fread("output/deseq2/ruakura/ruakura_analysis_sig_degs.csv")
not_concat_ru_ids <- not_concat_ru_degs$rn
concat_ru_ids <- rua_ordered_sig_res_group_table$`#gene_id`

vd <- venn.diagram(x = list("Concat Mapping Ru DEGs"=concat_ru_ids, "Regular Mapping Ru DEGs"=not_concat_ru_ids), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1, label=TRUE)
grid.newpage()
grid.draw(vd)
