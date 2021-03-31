library(data.table)
library(DESeq2)

sample_data <- fread("data/sample_table.csv", header=TRUE)

##########################
##Sample Parasitism Rate##
##########################

##Dunedin##
dunedin <- subset(sample_data, Weevil_Location=="Dunedin")
##percent parasitised
(sum(dunedin$Parasitism_status=="parasitised"))/(sum(dunedin$Weevil_Location=="Dunedin"))*100

dun_evasive <- subset(dunedin, Behaviour=="E")
(sum(dun_evasive$Parasitism_status=="parasitised"))/(sum(dun_evasive$Behaviour=="E"))*100

dun_nonev <- subset(dunedin, Behaviour=="N")
(sum(dun_nonev$Parasitism_status=="parasitised"))/(sum(dun_nonev$Behaviour=="N"))*100


##Ruakura##
ruakura <- subset(sample_data, Weevil_Location=="Ruakura")
##number parasitised
((sum(ruakura$Parasitism_status=="parasitised"))/(sum(ruakura$Weevil_Location=="Ruakura")))*100

rua_evasive <- subset(ruakura, Behaviour=="E")
(sum(rua_evasive$Parasitism_status=="parasitised"))/(sum(rua_evasive$Behaviour=="E"))

rua_nonev <- subset(ruakura, Behaviour=="N")
(sum(rua_nonev$Parasitism_status=="parasitised"))/(sum(rua_nonev$Behaviour=="N"))

##############################
##Plotting PCR Target Counts##
##############################

## ASW ##
asw_dds <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
asw_dds$group <- factor(paste(asw_dds$Parasitism_PCR,sep="_"))
##plot counts for PCR target gene - two ASW 'genes' COI in transcriptome
##TRINITY_DN63461_c12_g1
##TRINITY_DN7363_c1_g1
plotCounts(asw_dds, "ASW_TRINITY_DN7363_c1_g1", intgroup = c("group"))

## MH ##
mh_dds <- readRDS("output/deseq2/mh_dual/mh_dual_dds.rds")
mh_dds$group <- factor(paste(mh_dds$Parasitism_PCR,sep="_"))
##plot counts for PCR target gene
plotCounts(mh_dds, "MH_TRINITY_DN1053_c1_g1", intgroup = c("group"))
