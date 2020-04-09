library(data.table)
library(DESeq2)
library(fgsea)
library(ggplot2)
library(VennDiagram)

##LRT vs Wald
##LRT is good if you have multiple levels of a factor, so for example if you had three locations and you wanted to know if there was any interaction between location and parasitism. But in this case it's just a pairwise comparison, so you can test the additional effect of para vs. non_para to Ruakura vs. Lincoln.
##Wald is just a pairwise comparison, but if you use that you get interpretable L2FCs. We have to think about what the L2FC means for an interaction.
##wald FC is interpretable because it doesn't do fold-change shrinkage

##this ignores the fact some of the unpara asw were attacked though - do we trust the behaviour labels enough to alter them somewhat?

##########
##DESeq2##
##########

dds <- readRDS("output/deseq2/asw/asw_dds.rds")
dds_location_parasitism <- copy(dds)

dds_location_parasitism$Location <- factor(dds_location_parasitism$Weevil_Location)
dds_location_parasitism$Parasitism <- factor(dds_location_parasitism$Parasitism_status)
dds_location_parasitism$Sample_cleaned <- factor(dds_location_parasitism$Cleaned)

##Look for genes DE between locations, controlling for parasitism status
design(dds_location_parasitism) <- ~Sample_cleaned+Parasitism+Location+Parasitism:Location
dds_location_parasitism <- DESeq(dds_location_parasitism, test='Wald')
saveRDS(dds_location_parasitism, "output/deseq2/asw/WT_location-parasitism_int/asw_dds_WT_location-para-interaction.rds")
dds_location_parasitism <- readRDS("output/deseq2/asw/WT_location-parasitism_int/asw_dds_WT_location-para-interaction.rds")

##int term is last level - Parasitismundetected.LocationRuakura
resultsNames(dds_location_parasitism)


##so it's the additive effect of the second level of each factor
##it hasn't printed the base levels

##extract results -  lfc messed around by interaction terms - do not use LFC theshold with interactions
res_group <- results(dds_location_parasitism, name = c("Parasitismundetected.LocationRuakura"), alpha = 0.1)
##summarize results
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_res_group_table, "output/deseq2/asw/WT_location-parasitism_int/full_res.csv")
fwrite(ordered_sig_res_group_table, "output/deseq2/asw/WT_location-parasitism_int/sig_degs.csv")

##merge DEGs with trinotate annots
best_annot <- fread("data/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, best_annot, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw/WT_location-parasitism_int/sig_degs_annots.csv")

##Sub in any gene of interest to plot counts  
plotCounts(dds_location_parasitism, "ASW_TRINITY_DN12075_c0_g1", intgroup = c("Parasitism", "Location"), main="")

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
lp_ranks <- ordered_res_group_table[!is.na(stat), stat]
names(lp_ranks) <- ordered_res_group_table[!is.na(stat), rn]

lp_fgsea_res <- fgsea(pathways, ranks, nperm = 10000)
lp_sorted_fgsea_res <- fgsea_res[order(fgsea_res$padj)]
sum(sorted_fgsea_res$padj<0.05)


sig_fgsea_res <- subset(sorted_fgsea_res, padj < 0.1)
annot_sig_fgsea <- merge(sig_fgsea_res, go_annot_table, by.x="pathway", by.y="pathway", all.x=TRUE)


##to try plot all genes
topGenes <- rownames(ordered_sig_res_group_table)[order(ordered_sig_res_group_table$padj)][1:20]
z <- lapply(topGenes, function(x) plotCounts(dds_location_parasitism, x, c("Parasitism", "Location"),returnData = TRUE))
for(i in 1:20) z[[i]]$gene <- rep(topGenes[i], 79)
z <- do.call(rbind, z)
ggplot(z, aes(Location, count, colour = Parasitism)) + scale_y_log10() + geom_point(position = position_jitter(width = 0.1, height = 0), size = 3) + facet_wrap(~gene)



vd <- venn.diagram(x = list("WT location:parasitism"=wt_sig_names, "LRT location:attack"=sig_degs_names), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1)
grid.newpage()
grid.draw(vd)

