library(data.table)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(fgsea)

##LRT for unpara attacked vs unattacked
##do attacked but unpara asw differ from unattacked unpara asw?

##try without filtering for undetected so we can also include successful attack?##

##doesn't find anything of interest

dds <- readRDS("output/deseq2/asw/asw_dds.rds")
dds_LRT_location_attacked <- copy(dds)

##use to restrict to para asw
dds_LRT_location_attacked <- dds[,dds$Parasitism_status == "undetected"]


dds_LRT_location_attacked$Location <- factor(dds_LRT_location_attacked$Weevil_Location)
dds_LRT_location_attacked$Attacked <- factor(dds_LRT_location_attacked$Attacked)
dds_LRT_location_attacked$Sample_cleaned <- factor(dds_LRT_location_attacked$Cleaned)

##Look for genes DE between locations, controlling for parasitism status
design(dds_LRT_location_attacked) <- ~Sample_cleaned+Location+Attacked
dds_LRT_location_attacked <- DESeq(dds_LRT_location_attacked, test = "LRT", reduced=~Sample_cleaned+Location)
#saveRDS(dds_LRT_location_attacked, "output/deseq2/asw/LRT_attack/asw_dds_LRT_attack.rds")
#dds_LRT_location_attacked <- readRDS("output/deseq2/asw/LRT_attack/asw_dds_LRT_attack.rds")

resultsNames(dds_LRT_location_attacked)

##extract results -  lfc messed around by interaction terms - do not use LFC theshold with interactions
res_group <- results(dds_LRT_location_attacked, alpha = 0.1)
##summarize results
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
#fwrite(ordered_res_group_table, "output/deseq2/asw/LRT_attack/full_res.csv")
#fwrite(ordered_sig_res_group_table, "output/deseq2/asw/LRT_attack/sig_degs.csv")

plotCounts(dds_LRT_location_attacked, "ASW_TRINITY_DN13903_c0_g1", intgroup = c("Location", "Attacked"))

##merge DEGs with trinotate annots
best_annot <- fread("data/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, best_annot, by.x="rn", by.y="#gene_id", all.x=TRUE)
#fwrite(sig_annots, "output/deseq2/asw/LRT_attack/sig_degs_annots.csv")
