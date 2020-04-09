library(data.table)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(EnhancedVolcano)

dds_group <- readRDS("output/deseq2/dds_group.rds")

resultsNames(dds_group)

rua_res_group <- results(dds_group, contrast = c("group", "Ruakura_E", "Ruakura_N"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
rua_ordered_res_group <- rua_res_group[order(rua_res_group$padj),]
##Make data table and write to output
rua_ordered_res_group_table <- data.table(data.frame(rua_ordered_res_group), keep.rownames = TRUE)
rua_ordered_sig_res_group_table <- subset(rua_ordered_res_group_table, padj < 0.05)

##write full and sig res files
fwrite(rua_ordered_res_group_table, "output/deseq2/ruakura/res_group.csv")
fwrite(rua_ordered_sig_res_group_table, "output/deseq2/ruakura/ruakura_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)

##plot counts
plotCounts(dds_group, "TRINITY_DN32851_c0_g1", intgroup = c("group"), main = "TRINITY_DN41656_c0_g1 - RYamide receptor")
##volcano plot
EnhancedVolcano(rua_ordered_res_group_table, x="log2FoldChange", y="padj", lab="", pointSize = 3)

##read in annotated transcriptome
trinotate_most_sig_hits <- fread("data/asw_transcriptome/most_sig_transcript_blastx_hit_for_each_gene.csv")
setnames(rua_ordered_sig_res_group_table, old=c("rn"), new=c("#gene_id"))
##merge list of sig genes with annotations
rua_sig_w_annots <- merge(rua_ordered_sig_res_group_table, trinotate_most_sig_hits, by="#gene_id", all.x=TRUE)
##save file - no need to dedup anymore
fwrite(rua_sig_w_annots, "output/deseq2/ruakura/rua_sig_genes_with_annots.csv")

##sum of DEGs with no blastX annotation in transcriptome
sum(is.na(rua_sig_w_annots$sprot_Top_BLASTX_hit))
##list of DEGs with no blastX annotation
no_blastx_annot_degs <- rua_sig_w_annots[is.na(rua_sig_w_annots$sprot_Top_BLASTX_hit),]
##make list of degs with no blastx annot. - gene id but fasta is isoform ids - bbtools will filter for substring
list_degs_no_annot <- data.table(no_blastx_annot_degs$`#gene_id`)
##write list of degs with no blastx annot.
fwrite(list_degs_no_annot, "output/deseq2/ruakura/degs_with_no_annot.txt")

##blastx res for unann genes
blastx_unann_degs <- fread("output/deseq2/ruakura/unann/unann_degs_blastx.outfmt6")
setnames(blastx_unann_degs, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("#gene_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "unann_blastx_annotation"))
setorder(blastx_unann_degs, `#gene_id`, evalue, -bit_score)
blast_min_evalues <- blastx_unann_degs[,.SD[which.min(evalue)], by=`#gene_id`]
blast_min_evalues$`#gene_id` <- tstrsplit(blast_min_evalues$`#gene_id`, "_i", keep=c(1))
blast_min_evalues$unann_blastx_annotation <- tstrsplit(blast_min_evalues$unann_blastx_annotation, "<>", keep=c(1))
trinotate_blast_annots <- merge(rua_sig_w_annots, blast_min_evalues, by='#gene_id', all=TRUE)
fwrite(trinotate_blast_annots, "output/deseq2/ruakura/rua_degs_blast_trinotate.csv")

##DEGs with unchar or no annot - interproscan
unchar_or_hypo_annots <- dplyr::filter(trinotate_blast_annots, grepl('uncharacterized|hypothetical', unann_blastx_annotation))
##filter out genes with no manual annotation OR trinotate blastx annotation
no_manual_annot <- trinotate_blast_annots %>% filter(is.na(unann_blastx_annotation))
no_annot <- no_manual_annot[no_manual_annot$sprot_Top_BLASTX_hit == ".",]
##merge list of genes with no annot OR hypothetical/uncharacterised and save for interproscan
unchar_hypo_ids <- data.table(unchar_or_hypo_annots$transcript_id)
noannot_ids <- data.table(no_annot$transcript_id)
ids_for_interproscan <- merge(unchar_hypo_ids, noannot_ids, all = TRUE)
fwrite(ids_for_interproscan, "output/deseq2/ruakura/unann/interproscan_ids.txt")

