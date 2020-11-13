library(data.table)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(fgsea)
library(VennDiagram)

##do attacked but unpara asw differ from unattacked unpara asw?

dds <- readRDS("output/deseq2/asw/asw_dds.rds")
dds_location_attacked <- copy(dds)

##use to restrict to unpara
dds_location_attacked <- dds[,dds$Parasitism_status == "undetected"]
dds_location_attacked$Attacked <- factor(dds_location_attacked$Attacked)

##Look for genes DE between locations, controlling for parasitism status
design(dds_location_attacked) <- ~Attacked
dds_location_attacked <- DESeq(dds_location_attacked)
saveRDS(dds_location_attacked, "output/deseq2/asw/WT_attacked/asw_dds_WT_attacked_unpara.rds")
dds_location_attacked <- readRDS("output/deseq2/asw/WT_attacked/asw_dds_WT_attacked_unpara.rds")

resultsNames(dds_location_attacked)

##extract results -  lfc messed around by interaction terms - do not use LFC theshold with interactions
res_group <- results(dds_location_attacked, contrast = c("Attacked", "Y", "NO"), alpha = 0.1)
##summarize results
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_res_group_table, "output/deseq2/asw/WT_attacked/all_unpara_full_res.csv")
fwrite(ordered_sig_res_group_table, "output/deseq2/asw/WT_attacked/all_unpara_sig_degs.csv")

##use to plot with different sample groups
dds_location_attacked$Parasitism <- factor(dds_location_attacked$Parasitism_status)
dds_location_attacked$Location <- factor(dds_location_attacked$Weevil_Location)
##plot
plotCounts(dds_location_attacked, "ASW_TRINITY_DN4746_c0_g1", intgroup = c("Parasitism"))

##merge DEGs with trinotate annots
best_annot <- fread("data/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, best_annot, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw/WT_attacked/all_unpara_sig_degs_annots.csv")


###compare results from analyses
all_attacked_degs <- fread("output/deseq2/asw/WT_attacked/all_attacked_sig_degs_annots.csv")
unpara_attacked_degs <- fread("output/deseq2/asw/WT_attacked/sig_degs_annots.csv")
parasitism_degs <- fread("output/deseq2/asw/WT_parasitism/sig_degs_annots.csv")

vd <- venn.diagram(x = list("all attacked"=all_attacked_degs$rn, "unpara attacked"=unpara_attacked_degs$rn, "parasitism"=parasitism_degs$rn), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1,)
grid.newpage()
grid.draw(vd)

all_attacked_not_para_degs <- setdiff(all_attacked_degs$rn, parasitism_degs$rn)
all_attacked_not_para <- all_attacked_degs[all_attacked_degs$rn %in% all_attacked_not_para_degs,]
fwrite(all_attacked_not_para, "output/deseq2/asw/WT_parasitism/all_attacked_not_DE_parasitism_sig_annots.csv")

