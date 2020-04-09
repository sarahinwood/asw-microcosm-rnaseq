library(tximport)
library(data.table)
library(DESeq2)
library(VennDiagram)
library(dplyr)
library(EnhancedVolcano)

##need to make best hit per gene file for mh and edit this one's ids instead
trinotate_report <- fread("data/mh_edited_transcript_ids/trinotate_annotation_report.txt")

dds_concat_group_mh <- readRDS("output/concat_deseq2/mh_deseq2/dds_concat_group_mh.rds")

plotCounts(dds_concat_group_mh, "MH_TRINITY_DN1053_c1_g1", intgroup = c("group"))


#############
## Dunedin ##
#############
dun_viral_infection_res_group <- results(dds_concat_group_mh, contrast = c("group", "Dunedin_Yes", "Dunedin_No"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
dun_viral_infection_ordered_res_group <- dun_viral_infection_res_group[order(dun_viral_infection_res_group$padj),]
##Make data table and write to output
dun_viral_infection_ordered_res_group_table <- data.table(data.frame(dun_viral_infection_ordered_res_group), keep.rownames = TRUE)
dun_viral_infection_sig_res <- subset(dun_viral_infection_ordered_res_group_table, padj < 0.05)
dun_viral_infection_sig_res$rn <- tstrsplit(dun_viral_infection_sig_res$rn, "ASW_", keep=c(2))
dun_viral_inf_sig_annots <- merge(dun_viral_infection_sig_res, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(dun_viral_inf_sig_annots, "output/deseq2/virus_location/dun_viral_vs_non-viral_sig_annots.csv")
EnhancedVolcano(dun_viral_infection_sig_res, x="log2FoldChange", y="padj", lab="", pointSize = 3)

#############
## Ruakura ##
#############
rua_viral_infection_res_group <- results(dds_concat_group_mh, contrast = c("group", "Ruakura_Yes", "Ruakura_No"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
rua_viral_infection_ordered_res_group <- rua_viral_infection_res_group[order(rua_viral_infection_res_group$padj),]
##Make data table and write to output
rua_viral_infection_ordered_res_group_table <- data.table(data.frame(rua_viral_infection_ordered_res_group), keep.rownames = TRUE)
rua_viral_infection_sig_res <- subset(rua_viral_infection_ordered_res_group_table, padj < 0.05)
rua_viral_infection_sig_res$rn <- tstrsplit(rua_viral_infection_sig_res$rn, "ASW_", keep=c(2))
fwrite(rua_viral_infection_sig_res, "output/deseq2/virus_location/rua_viral_vs_non-viral_sig_degs.csv")
rua_viral_inf_sig_annots <- merge (rua_viral_infection_sig_res, trinotate_report, by.x="rn", by.y="#gene_id")
fwrite(rua_viral_inf_sig_annots, "output/deseq2/virus_location/rua_viral_vs_non-viral_sig_annots.csv")
EnhancedVolcano(rua_viral_infection_ordered_res_group_table, x="log2FoldChange", y="padj", lab="", transcriptPointSize = 3)
