library(tximport)
library(data.table)
library(DESeq2)
library(VennDiagram)
library(dplyr)
library(EnhancedVolcano)

trinotate_report <- fread("data/asw_transcriptome/most_sig_transcript_blastx_hit_for_each_gene.csv")
setnames(rua_ordered_sig_res_group_table, old=c("rn"), new=c("#gene_id"))

dds <- readRDS("output/deseq2/dds.rds")
##create dds object for group analysis
dds_group <- copy(dds)

##create groupings of weevil location and behavioural response to parasitoid
dds_group$group <- factor(paste(dds$Weevil_Location,dds$Behaviour,dds$Viral_expressed,sep="_"))

##add group to design
design(dds_group) <- ~group
##run deseq2 and generate results
dds_group <- DESeq(dds_group)
##save dds_group
saveRDS(dds_group, file = "output/deseq2/virus_location/dds_group_viral_location_behaviour.rds")

##create groupings of weevil location and viral expression
dds_viral <- copy(dds_group)
dds_viral$group <- factor(paste(dds$Weevil_Location,dds$Viral_expressed,sep="_"))
##add group to design
design(dds_viral) <- ~group
##run deseq2 and generate results
dds_viral <- DESeq(dds_viral)
##save dds_group
saveRDS(dds_viral, file = "output/deseq2/virus_location/dds_virus_location.rds")

#############
## Dunedin ##
#############
dun_viral_infection_res_group <- results(dds_viral, contrast = c("group", "Dunedin_Yes", "Dunedin_No"), lfcThreshold = 1, alpha = 0.1)
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
rua_viral_infection_res_group <- results(dds_viral, contrast = c("group", "Ruakura_Yes", "Ruakura_No"), lfcThreshold = 1, alpha = 0.1)
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

rua_ev_vs_nev <- fread("output/deseq2/ruakura/ruakura_analysis_sig_degs.csv")
rua_ev_vs_nev_ids <- rua_ev_vs_nev$rn

rua_viral_inf_deg_ids <- rua_viral_infection_sig_res$rn
dun_viral_inf_deg_ids <- dun_viral_infection_sig_res$rn

vd <- venn.diagram(x = list("Dunedin viral vs non-viral" = dun_viral_inf_deg_ids, "Ruakura viral vs non-viral"=rua_viral_inf_deg_ids), filename=NULL, cex = 1, cat.cex=1, lwd=1)
grid.newpage()
grid.draw(vd)

##rerun through with this design to see effect of viral gene expression without comparing locations separately??
##create groupings of weevil location and viral expression
dds_group$group <- factor(paste(dds$Viral_expressed))
##add group to design
design(dds_group) <- ~group
##run deseq2 and generate results
dds_group <- DESeq(dds_group)
##save dds_group
saveRDS(dds_group, file = "output/deseq2/virus_location/dds_group_viral.rds")

resultsNames(dds_viral)

viral_infection_res_group <- results(dds_group, contrast = c("group", "Yes", "No"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
viral_infection_ordered_res_group <- viral_infection_res_group[order(viral_infection_res_group$padj),]
##Make data table and write to output
viral_infection_ordered_res_group_table <- data.table(data.frame(viral_infection_ordered_res_group), keep.rownames = TRUE)
fwrite(viral_infection_ordered_res_group_table, "output/deseq2/virus_location/viral_expressed_all.csv")
viral_infection_ordered_sig_res_group_table <- subset(viral_infection_ordered_res_group_table, padj < 0.05)
viral_inf_sig_annots <- merge (viral_infection_ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id")
fwrite(viral_inf_sig_annots, "output/deseq2/virus_location/viral_vs_non-viral_sig_annots.csv")
EnhancedVolcano(viral_infection_ordered_res_group_table, x="log2FoldChange", y="padj", lab="", transcriptPointSize = 3)

plotCounts(dds_group, "TRINITY_DN4928_c0_g1", intgroup = c("group"), main="")

dedup_sig_annots <- fread("output/deseq2/virus_location/dedup_viral_vs_non-viral_sig_annots.csv")
sum(dedup_sig_annots$sprot_Top_BLASTX_hit==".")
##list of DEGs with no blastX annotation
no_blastx_annot_degs <- dedup_sig_annots[dedup_sig_annots$sprot_Top_BLASTX_hit == ".",]
##list of DEGs with no blastX OR blastP
no_blast_annot_degs <- no_blastx_annot_degs[no_blastx_annot_degs$sprot_Top_BLASTP_hit == ".",]
##make list of degs with no blast annot.
list_degs_no_annot <- data.table(no_blast_annot_degs$rn)
##write list of degs with no annot.
fwrite(list_degs_no_annot, "output/deseq2/virus_location/no_annot/degs_no_annot.txt")

#####read in extra BlastX results#####
##read in blast results
viral_no_annot_blast <- fread("output/deseq2/virus_location/no_annot/blastx.outfmt6")
##rename columns
setnames(viral_no_annot_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("#gene_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
setorder(viral_no_annot_blast, `#gene_id`, evalue, -bit_score)
fwrite(viral_no_annot_blast, "output/deseq2/virus_location/no_annot/blast_all_res.csv")

##extract result with lowest evalue for each peptide
blast_min_evalues <- viral_no_annot_blast[,.SD[which.min(evalue)], by=`#gene_id`]
blast_min_evalues$annotation <- tstrsplit(blast_min_evalues$annotation, "<>", keep=c(1))
blast_min_evalues$`#gene_id` <- tstrsplit(blast_min_evalues$`#gene_id`, "_i", keep=c(1))
fwrite(blast_min_evalues, "output/deseq2/virus_location/no_annot/blast_min_evalues.csv")

##merge blastx annotations for unann transcripts with transcriptome annotations
all_annots_degs <- merge(dedup_sig_annots, blast_min_evalues, by.x="rn", by.y="#gene_id", all = TRUE)
fwrite(all_annots_degs, "output/deseq2/virus_location/degs_trinotate_blastx_annots.csv")

##filter out genes in blastx annotation column that contain "uncharacterized" or "hypothetical"
unchar_or_hypo_annots <- dplyr::filter(all_annots_degs, grepl('uncharacterized|hypothetical', annotation))
##filter out genes with no manual annotation OR trinotate blastx annotation
no_manual_annot <- all_annots_degs %>% filter(is.na(annotation))
no_annot <- no_manual_annot[no_manual_annot$sprot_Top_BLASTX_hit == ".",]
##merge list of genes with no annot OR hypothetical/uncharacterised and save for interproscan
unchar_hypo_ids <- data.table(unchar_or_hypo_annots$transcript_id)
noannot_ids <- data.table(no_annot$transcript_id)
ids_for_interproscan <- merge(unchar_hypo_ids, noannot_ids, all = TRUE)
fwrite(ids_for_interproscan, "output/deseq2/virus_location/interproscan/interproscan_ids.txt")

##viral annots
trinotate_virus_annots <- dplyr::filter(all_annots_degs, grepl('virus', sprot_Top_BLASTX_hit))
manblast_virus_annots <- dplyr::filter(all_annots_degs, grepl('virus', annotation))
fwrite(virus_annots, "output/deseq2/virus_location/viral_annots.csv")
