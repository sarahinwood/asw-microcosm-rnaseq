library(data.table)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(fgsea)

##location-specific responses to attack status

dds <- readRDS("output/deseq2/asw/asw_dds.rds")
dds_location_attack <- copy(dds)

dds_location_attack$Location <- factor(dds_location_attack$Weevil_Location)
dds_location_attack$Attack <- factor(dds_location_attack$Attack_status)
dds_location_attack$Sample_cleaned <- factor(dds_location_attack$Cleaned)

##Look for genes DE between locations, controlling for attack status
design(dds_location_attack) <- ~Sample_cleaned+Attack+Location+Attack:Location
dds_location_attack <- DESeq(dds_location_attack, test='LRT', reduced=~Sample_cleaned+Attack+Location)
saveRDS(dds_location_attack, "output/deseq2/asw/LRT_location-attack-status_int/dds_int_location_attack.rds")

dds_location_attack <- readRDS("output/deseq2/asw/LRT_location-attack-status_int/dds_int_location_attack.rds")

##extract results
dds_res <- results(dds_location_attack, alpha=0.05)
summary(dds_res)
##results for all genes - for FGSEA analysis
full_res_dt <- data.table(data.frame(dds_res), keep.rownames=TRUE)
fwrite(full_res_dt, "output/deseq2/asw/LRT_location-attack-status_int/LRT_int_location_para_full_res.csv")

##make list of sig genes
sig_degs <- subset(dds_res, padj<0.05)
sig_degs_names <- row.names(sig_degs)
fwrite(data.table(sig_degs_names), "output/deseq2/asw/LRT_location-attack-status_int/sig_degs_names.csv")

##sig degs table
ordered_sig_degs <- sig_degs[order(sig_degs$padj),]
ordered_sig_degs_table <- data.table(data.frame(ordered_sig_degs), keep.rownames = TRUE)
ordered_sig_degs_table$fixed_ids <- tstrsplit(ordered_sig_degs_table$rn, "ASW_", keep=c(2))
fwrite(ordered_sig_degs_table, "output/deseq2/asw/LRT_location-attack-status_int/sig_degs.csv")
##merge with annots
trinotate_best_per_gene <- fread("data/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
ordered_sig_degs_annots <- merge(ordered_sig_degs_table, trinotate_best_per_gene, by.x="rn", by.y="#gene_id", all.x=TRUE, all.y=FALSE)
fwrite(ordered_sig_degs_annots, "output/deseq2/asw/LRT_location-attack-status_int/sig_degs_annots.csv")

##plot counts for genes of interest, sub in name
plotCounts(dds_location_attack, "ASW_TRINITY_DN3421_c0_g1", intgroup = c("Location", "Attack"))

dds_location_attack$loc_as <- factor(paste(dds_location_attack$Attack_status, dds_location_attack$Weevil_Location, sep="_"))
##plot expression pattern for gene
plot_gene <- plotCounts(dds_location_attack, "ASW_TRINITY_DN19594_c0_g1", 
                        intgroup = c("loc_as"), returnData = TRUE)
ggplot(plot_gene,
       aes(x = loc_as, y = count)) + 
  geom_boxplot() + geom_smooth(se = FALSE, method = "loess")+
  scale_y_log10() + ylab("Normalised Count")+
  xlab("Sample group")+
  ggtitle("")


##will plot all 20 top genes, but does so with set scale which makes some genes hard to see
topGenes <- (ordered_sig_degs_table$rn)[order(ordered_sig_degs_table$padj)][1:20]
z <- lapply(topGenes, function(x) plotCounts(dds_location_attack, x, c("Location", "Attack"),returnData = TRUE))
for(i in 1:20) z[[i]]$gene <- rep(topGenes[i], 79)
z <- do.call(rbind, z)
z <- within(z, joined <- paste(Location,Attack,sep= "_"))
z$joined <- z[paste(z$Location,z$Attack, sep="_")]
ggplot(z, aes(joined, count, colour=Attack)) + scale_y_log10() + geom_point(position = position_jitter(width = 0.1, height = 0), size = 3) + facet_wrap(~gene)

#########
##FGSEA##
#########
trinotate_report <- fread("data/asw_edited_transcript_ids/trinotate_annotation_report.txt", na.strings=".")
gene_ids <- trinotate_report[!is.na(gene_ontology_Pfam), unique(`#gene_id`)]
go_annot_list<-data.table(trinotate_report[,unique(unlist(strsplit(gene_ontology_Pfam, "`")))])
go_annot_table <- go_annot_list[,tstrsplit(V1, "^", fixed=TRUE)]
go_annot_table<-setnames(go_annot_table, old=c("V1", "V2", "V3"), new=c("pathway", "pathway_kind", "pathway_name"))

##function to extract GO terms from annotations in transcriptome (get all unique GO terms for each gene id) --> could look at other functional annot if I want to
EXTRACT_GO_TERMS <- function(x, trinotate_report){
  my_terms<-trinotate_report[`#gene_id`==x,unique(unlist(strsplit(gene_ontology_Pfam, "`")))]
  my_accessions<-unique(gsub("\\^.*", "", my_terms))
  my_accessions<-my_accessions[!is.na(my_accessions)]
  return(data.table(gene_id=x, accessions=my_accessions))
}

go_term_list <- lapply(gene_ids, EXTRACT_GO_TERMS, trinotate_report=trinotate_report)
go_term_table <- rbindlist(go_term_list)
term_to_gene <- go_term_table[,list(list(gene_id)), by=accessions]
pathways <- term_to_gene[,V1]
names(pathways) <- term_to_gene[,accessions]

##use stat column from deseq results to rank genes (can change if wanted)
setorder(full_res_dt, stat)
ranks <- full_res_dt[!is.na(stat), stat]
names(ranks) <- full_res_dt[!is.na(stat), rn]

fgsea_res <- fgsea(pathways, ranks, nperm = 10000)
sorted_fgsea_res <- fgsea_res[order(fgsea_res$padj)]
##44 enriched GO terms
sum(sorted_fgsea_res$padj<0.05)
sig_fgsea_res <- subset(sorted_fgsea_res, padj < 0.05)
annot_sig_fgsea <- merge(sig_fgsea_res, go_annot_table, by.x="pathway", by.y="pathway", all.x=TRUE)
fwrite(annot_sig_fgsea, "output/deseq2/asw/LRT_location-attack-status_int/sig_annot_fgsea_pfam.csv")
