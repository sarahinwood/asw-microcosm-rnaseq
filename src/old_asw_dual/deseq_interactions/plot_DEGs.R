library(data.table)
library(DESeq2)
library(tidyr)
library(ggplot2)

###plotcounts
dds <- readRDS("output/deseq2/asw/WT_location-parasitism_int/asw_dds_WT_location-para-interaction.rds")
dds$Location <- factor(dds$Weevil_Location)
dds$Attacked <- factor(dds$Attacked)
dds$Attack_status <- factor(dds$Attack_status)
dds$Parasitism <- factor(dds$Parasitism_status)
dds$location_para <- factor(paste(dds$Parasitism, dds$Location, sep="_"))
dds$loc_as <- factor(paste(dds$, dds$Location, sep="_"))
##normalized counts for plot
norm_counts <- counts(dds, normalized=TRUE)
norm_counts_dt <- data.table(norm_counts, keep.rownames = TRUE)

plot_degs <- fread("output/deseq2/asw/comp_interaction_res/DEGs_for_plots.txt", skip = c(1), header=FALSE)

##plot expression pattern for gene
plot_gene <- plotCounts(dds, "ASW_TRINITY_DN3390_c0_g1", 
                        intgroup = c("loc_as"), returnData = TRUE)
ggplot(plot_gene,
       aes(x = loc_as, y = count)) + 
  geom_point() + geom_smooth(se = FALSE, method = "loess")+
  scale_y_log10() + ylab("Normalised Count")+
  xlab("Sample group")+
  ggtitle("ASW_TRINITY_DN3390_c0_g1")


##
##look at location_attack status plot for this gene - not working currently
plotCounts(dds, "ASW_TRINITY_DN3390_c0_g1", 
           intgroup = c("loc_as"))




plot_gene_counts <- norm_counts_dt %>%
  filter(rn %in% plot_degs$V1)

plot_data <- plot_gene_counts %>%
  gather(colnames(plot_gene_counts)[2:9], key = "samplename", value = "normalized_counts")


ggplot(plot_data) +
  geom_point(aes(x = rn, y = normalized_counts)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

