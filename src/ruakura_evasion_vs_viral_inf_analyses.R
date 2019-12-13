library(data.table)
library(dplyr)
library(VennDiagram)
library(DESeq2)

rua_ev_vs_nev <- fread("output/deseq2/ruakura/ruakura_analysis_sig_degs.csv")
rua_ev_vs_nev_ids <- rua_ev_vs_nev$rn

rua_viral_inf <- fread("output/deseq2/virus_location/rua_viral_vs_non-viral_sig_degs.csv")
rua_viral_inf_deg_ids <- rua_viral_inf_deg_ids$rn

vd <- venn.diagram(x = list("Ruakura evasive
      vs non-evasive DEGs" = rua_ev_vs_nev_ids, "Ruakura viral vs non-viral"=rua_viral_inf_deg_ids), filename=NULL, cex = 1, cat.cex=1, lwd=1)
grid.newpage()
grid.draw(vd)

##find genes unique to ruakura evasion
evasion_specific_deg_ids <- setdiff(rua_ev_vs_nev_ids,rua_viral_inf_deg_ids)
evasion_specific_degs <- subset(rua_ev_vs_nev , (rn %in% evasion_specific_deg_ids))
##merge with annots
trinotate_report <- fread("data/asw_transcriptome/most_sig_transcript_blastx_hit_for_each_gene.csv")
rua_evasion_specific_annots <- merge(evasion_specific_degs, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(rua_evasion_specific_annots, "output/deseq2/virus_location/ru_behaviour_specific_degs.csv")


##do plotcount patterns of evasion specific degs match viral ones or do they look real?
dds_group <- readRDS("output/deseq2/dds_group.rds")
plotCounts(dds_group, "TRINITY_DN15411_c0_g2", intgroup = c("group"), main="T")

##overlap with nf exposed rnaseq
nf_exposed <- fread("output/nf_exposed/exposed_analysis_sig_degs.csv")
nf_exposed_ids <- nf_exposed$rn

vd2 <- venn.diagram(x = list("Ruakura Viral
      vs Non-viral DEGs" = rua_viral_inf_deg_ids, "Ruakura Ev vs N-Ev" = rua_ev_vs_nev_ids, "NF Exposed DEGs"=nf_exposed_ids), filename=NULL, cex = 1, cat.cex=1, lwd=1)
grid.newpage()
grid.draw(vd2)

##list 213 genes DE in all 3 analyses:
degs_all_behaviour_analyses <- intersect(rua_viral_inf_deg_ids, rua_ev_vs_nev_ids)
degs_all_analyses <- intersect(degs_all_behaviour_analyses, nf_exposed_ids)
degs_all_analyses_dt <- subset(rua_ev_vs_nev, (rn %in% degs_all_analyses))
degs_all_analyses_annots <- merge(degs_all_analyses_dt, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(degs_all_analyses_annots, "output/nf_exposed/degs_all_3_analyses_annots.csv")

                                