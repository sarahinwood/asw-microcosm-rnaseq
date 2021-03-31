library(data.table)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(EnhancedVolcano)
library(VennDiagram)

##Pairwise between dunedin and ruakura

dds <- readRDS("output/deseq2/asw/asw_dds.rds")

##create dds object for location analysis
dds_behaviour <- copy(dds)
dds_behaviour$behaviour <- factor(paste(dds$Behaviour))
dds_behaviour$group <- factor(paste(dds_behaviour$Behaviour, dds_behaviour$Weevil_Location, sep="_"))
dds_behaviour <- dds_behaviour[,dds_behaviour$Parasitism_status == "undetected"]

##add group to design
design(dds_behaviour) <- ~group
##run deseq2 and generate results
dds_behaviour <- DESeq(dds_behaviour)
##save dds_group
saveRDS(dds_behaviour, file = "output/deseq2/asw/WT_behaviour/asw_dds_behaviour.rds")

resultsNames(dds_behaviour)

##Make table of results for exposed vs control heads
res_group <- results(dds_behaviour, contrast = c("group", "E_Dunedin", "N_Dunedin"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
##write tables
fwrite(ordered_res_group_table, "output/deseq2/asw/WT_behaviour/full_res.csv")
fwrite(ordered_sig_res_group_table, "output/deseq2/asw/WT_behaviour/sig_degs.csv", col.names = TRUE, row.names = FALSE)
