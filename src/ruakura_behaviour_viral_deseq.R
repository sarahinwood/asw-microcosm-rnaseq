library(data.table)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(EnhancedVolcano)
library(VennDiagram)

dds_viral <- readRDS("output/deseq2/dds_viral.rds")
##read in annotated transcriptome
trinotate_most_sig_hits <- fread("data/asw_transcriptome/most_sig_transcript_blastx_hit_for_each_gene.csv")

#################################################################################
####DEG test between evasive and non-evasive for non-viral expressing samples####
#################################################################################

rua_no_res_group <- results(dds_viral, contrast = c("group", "Ruakura_E_No", "Ruakura_N_No"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
rua_no_ordered_res_group <- rua_no_res_group[order(rua_no_res_group$padj),]
##Make data table and write to output
rua_no_ordered_res_group_table <- data.table(data.frame(rua_no_ordered_res_group), keep.rownames = TRUE)
##nothing significant
rua_no_ordered_sig_res_group_table <- subset(rua_no_ordered_res_group_table, padj < 0.05)
##write full file
fwrite(rua_no_ordered_res_group_table, "output/deseq2/viral_ruakura/non-viral_res_group.csv")

#############################################################################
####DEG test between evasive and non-evasive for viral expressing samples####
#############################################################################

rua_yes_res_group <- results(dds_viral, contrast = c("group", "Ruakura_E_Yes", "Ruakura_N_Yes"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
rua_yes_ordered_res_group <- rua_yes_res_group[order(rua_yes_res_group$padj),]
##Make data table and write to output
rua_yes_ordered_res_group_table <- data.table(data.frame(rua_yes_ordered_res_group), keep.rownames = TRUE)
rua_yes_ordered_sig_res_group_table <- subset(rua_yes_ordered_res_group_table, padj < 0.05)

##write full and sig res files
fwrite(rua_yes_ordered_res_group_table, "output/deseq2/viral_ruakura/viral_res_group.csv")
fwrite(rua_yes_ordered_sig_res_group_table, "output/deseq2/viral_ruakura/viral_ruakura_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)

##merge with annotations
setnames(rua_yes_ordered_sig_res_group_table, old=c("rn"), new=c("#gene_id"))
##merge list of sig genes with annotations
rua_yes_sig_w_annots <- merge (rua_yes_ordered_sig_res_group_table, trinotate_most_sig_hits, by.x="#gene_id", by.y="#gene_id")
##save file - in excel edit duplicated gene ids (where one DEG had multiple annotations for each isoform)
fwrite(rua_yes_sig_w_annots, "output/deseq2/viral_ruakura/viral_rua_sig_genes_with_annots.csv")


plotCounts(dds_viral, "TRINITY_DN45187_c0_g1", intgroup = c("group"), main = "")


#########################################################################################
####compare overlap to analysis with all samples not controlling for viral expression####
#########################################################################################

no_viral_control_degs <- fread("output/deseq2/ruakura/ruakura_analysis_sig_degs.csv")
##make sig DEG lists
rua_viral_expressing <- rua_yes_ordered_sig_res_group_table$rn
rua_not_viral_expressing <- rua_no_ordered_sig_res_group_table$rn
rua_all <- no_viral_control_degs$rn
##make venn diagram
Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("Rua viral-expressing"=rua_viral_expressing, "Rua NOT viral-expressing"=rua_not_viral_expressing, "Rua all samples"=rua_all), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1, fill=Set1, label=TRUE)
grid.newpage()
grid.draw(vd)

