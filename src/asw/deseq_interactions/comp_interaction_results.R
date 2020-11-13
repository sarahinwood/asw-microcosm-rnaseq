library(data.table)
library(VennDiagram)

##overlap with location:parasitism analysis?
location_para_degs <- fread("output/deseq2/asw/WT_location-parasitism_int/sig_degs.csv")
location_attack_status_degs <- fread("output/deseq2/asw/LRT_location-attack-status_int/sig_degs.csv")
location_attacked_degs <- fread("output/deseq2/asw/WT_location-attacked_int/sig_degs.csv")
##make table of all genes and their res (if DE) in each comparison
first_merge <- merge(location_attacked_degs, location_attack_status_degs, by="rn", all=TRUE)
all_int_degs <- merge(first_merge, location_para_degs, by="rn", all=TRUE)
##merge with annots
trinotate_longest <- fread("data/asw_edited_transcript_ids/trinotate_longest_isoform.csv")
all_int_degs_annots <- merge(all_int_degs, trinotate_longest, by.x="rn", by.y="#gene_id")
setorder(all_int_degs_annots, padj.x, padj.y, padj, na.last=TRUE)
fwrite(all_int_degs_annots, "output/deseq2/asw/comp_interaction_res/all_interaction_degs.csv")

vd <- venn.diagram(x = list("Location:Parasitism"=location_para_degs$rn, "Location:Attack_status"=location_attack_status_degs$rn, "Location:Attacked"=location_attacked_degs$rn ), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1,)
grid.newpage()
grid.draw(vd)

##get lists of genes in each region
#intersect
#setdiff
#union
##2 genes shared by all
para_attacked <- intersect(location_para_degs$rn, location_attacked_degs$rn)
##16 shared between
para_attackstatus <- intersect(location_para_degs$rn, location_attack_status_degs$rn)

##location attacked specific
location_attacked <- data.table(setdiff(location_attacked_degs$rn, location_attack_status_degs$rn))
location_attacked_annots <- merge(location_attacked, trinotate_longest, by.x="V1", by.y="#gene_id")


##location_parasitism specific
location_parasitism <- data.table(setdiff(location_para_degs$rn, location_attack_status_degs$rn))
location_parasitism_annots <- merge(location_parasitism, trinotate_longest, by.x="V1", by.y="#gene_id")


##location_attack_status specific
location_attack_status <- setdiff(location_attack_status_degs$rn, location_para_degs$rn)
location_attack_status <- data.table(setdiff(location_attack_status, location_attacked_degs$rn))
location_attack_status_annots <- merge(location_attack_status, trinotate_longest, by.x="V1", by.y="#gene_id")

##misses ASW_TRINITY_DN5454_c0_g1 and ASW_TRINITY_DN12072_c0_g1 - why?
shared <- intersect(ordered_sig_degs_annots$rn, location_para_degs$rn)
