library(tximport)
library(data.table)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(viridis)
library(tidyverse)

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
##factors and design
asw_dds$para <- factor(paste(asw_dds$Parasitism_status))
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))
asw_dds$pc1_sign <- factor(paste(asw_dds$PC1_sign))
##relevel factors
asw_dds$para <- relevel(asw_dds$para, ref="undetected")
asw_dds$location <- relevel(asw_dds$location, ref="Dunedin")

asw_dds_para <- copy(asw_dds)
design(asw_dds_para) <- ~pc1_sign+location+para+location:para
##run deseq2
asw_dds_para <- DESeq(asw_dds_para)
saveRDS(asw_dds_para, "output/deseq2/asw_dual/INT_WT_parasitism-location/parasitism-location-int_WT.rds")

asw_dds_para <- readRDS("output/deseq2/asw_dual/INT_WT_parasitism-location/parasitism-location-int_WT.rds")
##results
res_group <- results(asw_dds_para, lfcThreshold = 1, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/asw_dual/INT_WT_parasitism-location/sig_degs.csv")
fwrite(ordered_res_group_table, "output/deseq2/asw_dual/INT_WT_parasitism-location/res_group.csv")

asw_dds_para$group <- factor(paste(asw_dds_para$Parasitism_status, asw_dds_para$location, sep=" "))
##plot
plot <- plotCounts(asw_dds_para, "ASW_TRINITY_DN30769_c1_g2", intgroup=("group"), returnData = TRUE)
ggplot(plot, aes(x=group, y=count))+
  ylab("Normalised counts")+
  xlab("")+
  geom_point(size=3, alpha=0.5, show.legend = FALSE, colour="#440154FF")+
  theme_bw()

trinotate_report <- fread("data/asw-mh-combined-transcriptome/output/asw_edited_transcript_ids/trinotate_longest_isoform.csv", na.strings = ".")
sig_annots <- merge(ordered_sig_res_group_table, trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)
fwrite(sig_annots, "output/deseq2/asw_dual/INT_WT_parasitism-location/sig_annots.csv")

##########
## plot ##
##########

##get gene counts
counts_table <- data.table(counts(asw_dds_para, normalized=TRUE), keep.rownames = TRUE)
annot_counts <- filter(counts_table, rn %in% sig_annots$rn)
##melt for plotting
plot_annots_counts <- annot_counts %>% gather(colnames(annot_counts)[2:80], key="sample_name", value="normalized_counts")
##sample group information
sample_to_group <- data.frame(colData(asw_dds_para)[,c(1,9,12)])
sample_to_group$group <- paste(sample_to_group$Weevil_Location, sample_to_group$Parasitism_status, sep="_")
plotting_counts <- inner_join(plot_annots_counts, sample_to_group)
##plot all annot DEGs using ggplot2
ggplot(plotting_counts) +
  geom_point(aes(x = group, y = normalized_counts, colour=group), alpha=0.5) +
  labs(colour="Group", y="Normalized counts", x="")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~rn, scales="free")

