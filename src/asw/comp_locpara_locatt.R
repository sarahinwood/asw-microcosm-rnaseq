library(data.table)
library(VennDiagram)

##overlap with location:parasitism analysis?
location_para_degs <- fread("output/deseq2/asw/WT_location-parasitism_int/sig_degs_annots.csv")
location_attack_degs <- fread("output/deseq2/asw/LRT_int_location_attack/sig_degs_annots.csv")


location_attack <- setdiff(ordered_sig_degs_annots$rn, location_para_degs$rn)
loc_att_sp <- location_attack_degs[location_attack_degs$rn %in% location_attack,]


location_parasitism <- setdiff(location_para_degs$rn, ordered_sig_degs_annots$rn)
loc_para_sp <- location_para_degs[location_para_degs$rn %in% location_parasitism,]

##misses ASW_TRINITY_DN5454_c0_g1 and ASW_TRINITY_DN12072_c0_g1 - why?

shared <- intersect(ordered_sig_degs_annots$rn, location_para_degs$rn)

vd <- venn.diagram(x = list("Location:Parasitism"=location_para_degs$rn, "Location:Attack"=location_attack_degs$rn), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1,)
grid.newpage()
grid.draw(vd)


plotCounts(dds_location_attack, "ASW_TRINITY_DN16220_c0_g1", intgroup = c("Location", "Attack"))

