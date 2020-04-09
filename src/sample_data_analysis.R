library("data.table")

sample_data <- fread("data/sample_key.csv", header=TRUE)

##########################
##Sample Parasitism Rate##
##########################

##Dunedin##
dunedin <- subset(sample_data, Weevil_Location=="Dunedin")
##number parasitised
(sum(dunedin$Parasitism_status=="parasitised"))/(length(dunedin$Sample_name))

dun_evasive <- subset(dunedin, Behaviour=="E")
(sum(dun_evasive$Parasitism_status=="parasitised"))/(length(dunedin$Sample_name))

dun_nonev <- subset(dunedin, Behaviour=="N")
(sum(dun_nonev$Parasitism_status=="parasitised"))/(length(dunedin$Sample_name))


##Ruakura##
ruakura <- subset(sample_data, Weevil_Location=="Ruakura")
##number parasitised
(sum(ruakura$Parasitism_status=="parasitised"))/(length(ruakura$Sample_name))

rua_evasive <- subset(ruakura, Behaviour=="E")
(sum(rua_evasive$Parasitism_status=="parasitised"))/(length(dunedin$Sample_name))

rua_nonev <- subset(ruakura, Behaviour=="N")
(sum(rua_nonev$Parasitism_status=="parasitised"))/(length(dunedin$Sample_name))

##############################
##Plotting PCR Target Counts##
##############################

##ASW##
asw_dds <- readRDS("output/deseq2/asw/asw_dds.rds")
asw_dds$group <- factor(paste(asw_dds$Parasitism_PCR,sep="_"))
##plot counts for PCR target gene
plotCounts(asw_dds, "ASW_TRINITY_DN1031_c1_g2", intgroup = c("group"))

##MH##

mh_dds <- readRDS("output/deseq2/mh/mh_dds.rds")
mh_dds$group <- factor(paste(mh_dds$Parasitism_PCR,sep="_"))
##plot counts for PCR target gene
plotCounts(mh_dds, "ASW_TRINITY_DN1031_c1_g2", intgroup = c("group"))
