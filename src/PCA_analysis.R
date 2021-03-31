library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)

#########
## ASW ##
#########

asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
##get factors ready for plotting
asw_dds$para_pcr <- factor(paste(asw_dds$Parasitism_PCR))
asw_dds$RQN <- factor(paste(asw_dds$RQN))
asw_dds$concentration <- factor(paste(asw_dds$conc))
asw_dds$location <- factor(paste(asw_dds$Weevil_Location))

##VST transformation
asw_vst <- varianceStabilizingTransformation(asw_dds, blind=TRUE)
##plot PCA with first 2 dimensions to investigate sample clustering
pca_plot <- plotPCA(asw_vst, intgroup=c("para_pcr"), returnData=TRUE)
percentVar <- round(100 * attr(pca_plot, "percentVar"))

##PCA plot (save with dim.s 3.00 x 8.00)
ggplot(pca_plot, aes(x=PC1, y=PC2, color=para_pcr))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="Parastism multiplex\n RT-PCR result")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()

##PCA with concentration and RQN
asw_pca_table <- data.table(plotPCA(asw_vst, intgroup=c("RQN"), returnData=TRUE))
##conc
asw_pca_table$concentration <- sapply(asw_pca_table$concentration, as.character)
asw_pca_table$concentration <- sapply(asw_pca_table$concentration, as.numeric)
##RQN
asw_pca_table$RQN <- sapply(asw_pca_table$RQN, as.character)
asw_pca_table$RQN <- sapply(asw_pca_table$RQN, as.numeric)

##PCA plot (save with dim.s 3.00 x 8.00)
##RNA\nconcentration\n(ng/uL) for concentration
ggplot(asw_pca_table, aes(x=PC1, y=PC2, color=RQN))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis()+
  labs(colour="RQN")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()

########
## Mh ##
########

mh_dds <- readRDS("output/deseq2/mh_dual/mh_dual_dds.rds")
##get factors ready for plotting
mh_dds$para_pcr <- factor(paste(mh_dds$Parasitism_PCR))
mh_dds$RQN <- factor(paste(mh_dds$RQN))
mh_dds$concentration <- factor(paste(mh_dds$conc))
##VST transformation
mh_vst <- varianceStabilizingTransformation(mh_dds, blind=TRUE)

##PCA plotting
mh_pca_plot <- plotPCA(mh_vst, intgroup=c("para_pcr"), returnData=TRUE)
percentVar <- round(100 * attr(mh_pca_plot, "percentVar"))
##PCA plot (save with dim.s 3.00 x 8.00)
ggplot(mh_pca_plot, aes(x=PC1, y=PC2, color=para_pcr))+
  geom_point(size=3, alpha=0.7)+
  scale_color_viridis(discrete=TRUE)+
  labs(colour="Parastism multiplex\n RT-PCR result")+
  xlab(paste("PC1:", percentVar[1], "% variance")) + 
  ylab(paste("PC2:", percentVar[2], "% variance")) + 
  coord_fixed()+
  theme_bw()
