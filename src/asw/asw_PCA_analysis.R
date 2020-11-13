library("DESeq2")
library("data.table")
library("ggplot2")

sample_data <- fread("data/sample_key.csv")
##read in dds saved in previous script
dds <- readRDS("output/deseq2/asw/asw_dds.rds")

##PCA plot - location_behaviour (first must log transform using vst (must set blind = true))
##makes expression comparable between samples?
vst <- varianceStabilizingTransformation(dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
##can change intgroup to anoy factor/s in sample key to change point colours
plotPCA(vst, intgroup=c("Weevil_Location"))

PCA_dt <- plotPCA(vst, intgroup=c("Parasitism_PCR"), returnData=TRUE)
fwrite(PCA_dt, "output/deseq2/asw/PCA/PCA_dt.csv")

####################################################################################
###not yet edited below here###

##viral expression
vst_viral <- varianceStabilizingTransformation(dds_viral, blind=TRUE)
plotPCA(vst_viral, intgroup=c("Viral_expressed"))
viral_PCA <- plotPCA(vst_viral, intgroup=c("Viral_expressed"), returnData = TRUE)
fwrite(viral_PCA, "output/deseq2/viral_PCA_table.csv")

##make table of PCA and viral counts data
viral_counts <- fread("output/deseq2/viral_counts_matrix.csv")
melted_viral_counts <- melt(viral_counts)
viral_counts_PCA <- merge(viral_PCA, melted_viral_counts, by.x="name", by.y="variable")
fwrite(viral_counts_PCA, "output/deseq2/viral_counts_PCA.csv")
viral_counts_PCA_sample_data <- merge(sample_data, viral_counts_PCA, by.x="Sample_name", by.y="name")
##use this file to check that viral expression labels for samples are correct
fwrite(viral_counts_PCA_sample_data, "output/deseq2/sample_data_viral_counts_PCA.csv")

#######edit for behavioural samples
##design formula??
##      ~Weevil_Location+Behaviour

dds_viral <- readRDS("output/deseq2/dds_viral.rds")

##pull out results from dds
rua_res_group <- results(dds, contrast = c("group", "Ruakura_E", "Ruakura_N"), lfcThreshold = 1, alpha = 0.1)
##only keep genes where padj is not NA
kept_genes <- rownames(subset(rua_res_group, !is.na(padj)))
##create matrix of vst values for only genes where padj didn't = NA
vst_asssay<- assay(vst)[kept_genes,]
##perform PCA on vst data matrix
pc <- prcomp(t(vst_asssay), center = TRUE, scale = TRUE)
##generate data table of results
pc_wide <- data.table(pc$x, keep.rownames = TRUE)
pc_pd <- melt(pc_wide)
fwrite(pc_pd, "output/deseq2/PCA/vst_pca_plot_data.csv")

pc_pd_sample_data <- merge(pc_pd, sample_data, by.x="rn", by.y="Sample_name")

##plot pcs
ggplot(pc_pd_sample_data, aes(x=rn, y=value, colour=Viral_expressed, shape=Weevil_Location))+
  facet_wrap(~variable)+geom_point()+theme(axis.text.x=element_text(angle = 90))

colData(dds)

###########effect of expressing viral genes on PCA###########

dds <- readRDS("output/deseq2/dds.rds")
vst <- varianceStabilizingTransformation(dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
plotPCA(vst, intgroup=c("Behaviour"))
