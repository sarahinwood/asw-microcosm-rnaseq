library(data.table)
library(DESeq2)
library(fgsea)
library(ggplot2)
library(VennDiagram)

##LRT vs Wald
##LRT is good if you have multiple levels of a factor, so for example if you had three locations and you wanted to know if there was any interaction between location and parasitism. But in this case it's just a pairwise comparison, so you can test the additional effect of para vs. non_para to Ruakura vs. Lincoln.
##Wald is just a pairwise comparison, but if you use that you get interpretable L2FCs. We have to think about what the L2FC means for an interaction.
##wald FC is interpretable because it doesn't do fold-change shrinkage

##########
##DESeq2##
##########

dds <- readRDS("output/deseq2/asw/asw_dds.rds")
dds_location_attacked <- copy(dds)

dds_location_attacked$Location <- factor(dds_location_attacked$Weevil_Location)
dds_location_attacked$Attack <- factor(dds_location_attacked$Attacked)
dds_location_attacked$Sample_cleaned <- factor(dds_location_attacked$Cleaned)
dds_location_attacked$Parasitism <- factor(dds_location_attacked$Parasitism_status)

##Look for genes DE between locations, controlling for parasitism status
design(dds_location_attacked) <- ~Sample_cleaned+Parasitism+Attack+Location+Attack:Location
dds_location_attacked <- DESeq(dds_location_attacked, test='Wald')
saveRDS(dds_location_attacked, "output/deseq2/asw/WT_location-attacked_int/asw_dds_WT_location-attacked_interaction.rds")
dds_location_attacked <- readRDS("output/deseq2/asw/WT_location-attacked_int/asw_dds_WT_location-attacked_interaction.rds")

##int term is last level - Parasitismundetected.LocationRuakura
resultsNames(dds_location_attacked)


##so it's the additive effect of the second level of each factor
##it hasn't printed the base levels

##extract results -  lfc messed around by interaction terms - do not use LFC theshold with interactions
res_group <- results(dds_location_attacked, name = c("AttackY.LocationRuakura"), alpha = 0.1)
##summarize results
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_res_group_table, "output/deseq2/asw/WT_location-attacked_int/full_res.csv")
fwrite(ordered_sig_res_group_table, "output/deseq2/asw/WT_location-attacked_int/sig_degs.csv")

##write sig gene nme list for clustering
fwrite(list(ordered_sig_res_group_table$rn), "output/deseq2/asw/WT_location-attacked_int/sig_gene_names.csv")

##merge DEGs with trinotate annots
best_annot <- fread("data/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, best_annot, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw/WT_location-attacked_int/sig_degs_annots.csv")

##Sub in any gene of interest to plot counts  
plotCounts(dds_location_attacked, "ASW_TRINITY_DN3390_c0_g1", intgroup = c("Attack", "Location"), main="")

dds_location_attacked$loc_att <- factor(paste(dds_location_attacked$Weevil_Location, dds_location_attacked$Attacked, sep="_"))
##plot expression pattern for gene
plot_gene <- plotCounts(dds_location_attacked, "ASW_TRINITY_DN2983_c1_g1", 
                        intgroup = c("loc_att"), returnData = TRUE)
ggplot(plot_gene,
       aes(x = loc_att, y = count)) + 
  geom_point() + geom_smooth(se = FALSE, method = "loess")+
  scale_y_log10() + ylab("Normalised Count")+
  xlab("Sample group")+
  ggtitle("ASW_TRINITY_DN2983_c1_g1 - 2-(3-amino-3-carboxypropyl)histidine synthase subunit 2")


#########
##FGSEA##
#########
trinotate_report <- fread("data/asw_edited_transcript_ids/trinotate_annotation_report.txt", na.strings = ".")
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
setorder(ordered_res_group_table, stat)
ranks <- ordered_res_group_table[!is.na(stat), stat]
names(ranks) <- ordered_res_group_table[!is.na(stat), rn]

fgsea_res <- fgsea(pathways, ranks, nperm = 10000)
sorted_fgsea_res <- fgsea_res[order(fgsea_res$padj)]
sum(sorted_fgsea_res$padj<0.05)


sig_fgsea_res <- subset(sorted_fgsea_res, padj < 0.1)
annot_sig_fgsea <- merge(sig_fgsea_res, go_annot_table, by.x="pathway", by.y="pathway", all.x=TRUE)
fwrite(annot_sig_fgsea, "output/deseq2/asw/WT_location-attacked_int/fgsea_sig_annot.csv")

