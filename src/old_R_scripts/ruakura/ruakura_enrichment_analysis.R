library(data.table)
library(fgsea)
library(ggplot2)
library(viridis)

##can search for blast and replace with pfam and v.v.

trinotate_report <- fread("data/asw_transcriptome/trinotate_annotation_report.txt", na.strings = ".")
gene_ids <- trinotate_report[!is.na(gene_ontology_blast), unique(`#gene_id`)]
res_group <- fread("output/deseq2/ruakura/res_group.csv")

go_annot_list<-data.table(trinotate_report[,unique(unlist(strsplit(gene_ontology_blast, "`")))])
go_annot_table <- go_annot_list[,tstrsplit(V1, "^", fixed=TRUE)]
go_annot_table<-setnames(go_annot_table, old=c("V1", "V2", "V3"), new=c("pathway", "pathway_kind", "pathway_name"))

##function to extract GO terms from annotations in transcriptome (get all unique GO terms for each gene id) --> could look at other functional annot if I want to
EXTRACT_GO_TERMS <- function(x, trinotate_report){
  my_terms<-trinotate_report[`#gene_id`==x,unique(unlist(strsplit(gene_ontology_blast, "`")))]
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
setorder(res_group, stat)
ranks <- res_group[!is.na(stat), stat]
names(ranks) <- res_group[!is.na(stat), rn]

fgsea_res <- fgsea(pathways, ranks, nperm = 10000)
sorted_fgsea_res <- fgsea_res[order(fgsea_res$padj)]

##function to calculate length of leading edge for each go term
go_terms <- sorted_fgsea_res$pathway
LEADING_EDGE_LENGTH <- function(x, sorted_fgsea_res){
  term_res <- fgsea_res[fgsea_res$pathway == x,]
  term_leading_edge <- data.frame(term_res$leadingEdge)
  leading_edge_length <- nrow(term_leading_edge)
  return(data.table(pathway=x, leading_edge_length=leading_edge_length))
}
LE_length <- lapply(go_terms, LEADING_EDGE_LENGTH, sorted_fgsea_res=sorted_fgsea_res)
LE_length_dt <- rbindlist(LE_length)
fgsea_res_full <- merge(fgsea_res, LE_length_dt, by="pathway")

sum(fgsea_res_full$padj<0.05)
fwrite(fgsea_res_full, "output/fgsea/ruakura/blast/ruakura_e_vs_n_fgsea_GOtermblast_deseqstat_res.csv")

##subset into only sig terms and merge w/annotations
sig_fgsea_res <- subset(fgsea_res_full, padj < 0.1)
annot_sig_fgsea <- merge(sig_fgsea_res, go_annot_table, by.x="pathway", by.y="pathway", all.x=TRUE)
fwrite(annot_sig_fgsea, "output/fgsea/ruakura/blast/ruakura_e_vs_n_sig_annot_fgsea_blast.csv")
##split into 3 tables --> biological process, cellular component and molecular function
bp_res <- annot_sig_fgsea[annot_sig_fgsea$`pathway_kind`=="biological_process"]
cc_res <- annot_sig_fgsea[annot_sig_fgsea$`pathway_kind`=="cellular_component"]
mf_res <- annot_sig_fgsea[annot_sig_fgsea$`pathway_kind`=="molecular_function"]

##plot normalised enrichment for GO terms where padj<0.1 (but indicate if padj<0.05) - can change to only bp, cc or mf
ggplot(bp_res, aes(reorder(pathway_name, NES), NES)) +
  geom_text(aes(y= NES + 0.3*sign(NES), label=round(padj, digits=3)), vjust=0, hjust=0) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Biological Process GO Pathway", y="FGSEA Normalized Enrichment Score") + 
  theme_minimal()+
  theme(axis.text.y = element_text(size=15))

##fgsea plot with dot position  = NES, colour = padj, size = leading edge size
## - probably better when not all negatively enriched
##geom-point twice to get plot reordered but also have point on top of grey line
##can use this as colour for line to easily show terms padj<0.05 - colour=ifelse((bp_res$padj<0.05), "black", "grey")
##biological process
ggplot(bp_res, aes(reorder(pathway_name, NES), NES)) + 
  geom_point(aes(size=bp_res$leading_edge_length, col=padj))  +
  geom_segment(aes(y = 0, x = pathway_name, yend = NES, xend = pathway_name), colour="darkgrey") +
  geom_point(aes(size=bp_res$leading_edge_length, col=padj))  +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_viridis() +
  theme_minimal()+
  coord_flip()+
  labs(x="Biological Process GO Pathway", y="FGSEA Normalized Enrichment Score",
  size="No. genes in
leading edge", col="Padj")

##molecular function
ggplot(mf_res, aes(reorder(pathway_name, NES), NES)) + 
  geom_point(aes(size=mf_res$leading_edge_length, col=padj))  +
  geom_segment(aes(y = 0, x = pathway_name, yend = NES, xend = pathway_name), colour="darkgrey") +
  geom_point(aes(size=mf_res$leading_edge_length, col=padj))  +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_viridis() +
  theme_minimal()+
  coord_flip()+
  labs(x="Molecular Function GO Pathway", y="FGSEA Normalized Enrichment Score",
       size="No. genes in
leading edge", col="Padj")

##cell comp.
ggplot(cc_res, aes(reorder(pathway_name, NES), NES)) + 
  geom_point(aes(size=cc_res$leading_edge_length, col=padj))  +
  geom_segment(aes(y = 0, x = pathway_name, yend = NES, xend = pathway_name), colour="darkgrey") +
  geom_point(aes(size=cc_res$leading_edge_length, col=padj))  +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_viridis() +
  theme_minimal()+
  coord_flip()+
  labs(x="Cellular Component GO Pathway", y="FGSEA Normalized Enrichment Score",
       size="No. genes in
leading edge", col="Padj")

##plot enrichment of GO term
plotEnrichment(pathways[["GO:0098586"]], ranks) + labs(title="cellular response to virus")

##GSEA table plot
topPathways <- fgsea_res[head(order(padj), n=15)][order(NES), pathway]
plotGseaTable(pathways[topPathways], ranks, fgsea_res, gseaParam = 0.5)

####find CORE members that contribute to ES score (present in list before running sum reaches max.dev. from 0)
core_mem_res <- fgsea_res[fgsea_res$pathway == "GO:0004983",]
term_leading_edge <- data.frame(core_mem_res$leadingEdge)
setnames(term_leading_edge, old=c("c..TRINITY_DN41656_c0_g1....TRINITY_DN9760_c0_g1.."), new=c("gene_id"))
term_leading_annots <- merge(term_leading_edge, trinotate_report, by.x="gene_id", by.y="#gene_id")
fwrite(term_leading_annots, "output/fgsea/ruakura/neuropeptideY_leading_edge_annots.csv")
##plot enrichment of GO term
plotEnrichment(pathways[["GO:0003723"]], ranks) + labs(title="neuropeptide Y")

