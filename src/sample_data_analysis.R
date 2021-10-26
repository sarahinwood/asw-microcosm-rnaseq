library(data.table)
library(DESeq2)

sample_data <- fread("data/sample_table.csv", header=TRUE)

##############################
##Plotting PCR Target Counts##
##############################

## ASW ##
asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
asw_dds$group <- factor(paste(asw_dds$Parasitism_PCR))
##plot counts for PCR target gene - two ASW 'genes' COI in transcriptome
##ASW_TRINITY_DN63461_c12_g1
##ASW_TRINITY_DN7363_c1_g1
plotCounts(asw_dds, "ASW_TRINITY_DN7363_c1_g1", intgroup = c("group"))

## MH ##
mh_dds <- readRDS("output/deseq2/mh_dual/mh_dual_dds.rds")
mh_dds$group <- factor(paste(mh_dds$Parasitism_PCR))
##plot counts for PCR target gene
##MH_TRINITY_DN481_c1_g1
##MH_TRINITY_DN17214_c0_g2 - low counts
Mh_PCR_counts <- plotCounts(mh_dds, "MH_TRINITY_DN481_c1_g1", intgroup = c("group"), returnData = TRUE)

ggplot(Mh_PCR_counts, aes(x=group, y=count, colour=group))+
  geom_point(size=3, alpha=0.7)+
  labs(y="Normalized counts", x="Parasitism PCR result")+
  scale_colour_viridis(discrete=TRUE)+
  scale_y_continuous(trans="log10")+
  theme_bw()+
  theme(legend.position="none")

##########################
##Sample Parasitism Rate##
##########################

##Dunedin##
dunedin <- subset(sample_data, Weevil_Location=="Dunedin")
##percent Dunedin parasitised
(sum(dunedin$Parasitism_status=="parasitized"))/(sum(dunedin$Weevil_Location=="Dunedin"))*100
##percent evasive parasitised
dun_evasive <- subset(dunedin, Behaviour=="E")
(sum(dun_evasive$Parasitism_status=="parasitized"))/(sum(dun_evasive$Behaviour=="E"))*100
##percent non-evasive parasitised
dun_nonev <- subset(dunedin, Behaviour=="N")
(sum(dun_nonev$Parasitism_status=="parasitized"))/(sum(dun_nonev$Behaviour=="N"))*100

##Ruakura##
ruakura <- subset(sample_data, Weevil_Location=="Ruakura")
##number parasitised
((sum(ruakura$Parasitism_status=="parasitized"))/(sum(ruakura$Weevil_Location=="Ruakura")))*100
##percent evasive parasitised
rua_evasive <- subset(ruakura, Behaviour=="E")
(sum(rua_evasive$Parasitism_status=="parasitized"))/(sum(rua_evasive$Behaviour=="E"))
##percent non-evasive parasitised
rua_nonev <- subset(ruakura, Behaviour=="N")
(sum(rua_nonev$Parasitism_status=="parasitized"))/(sum(rua_nonev$Behaviour=="N"))


