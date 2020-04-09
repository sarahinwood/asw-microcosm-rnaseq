library("data.table")
library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("EnhancedVolcano")

##combine location + behaviour
##not v. helpful due to para statuses of sample groups

##what if we remove all para samples?

dds <- readRDS("output/deseq2/asw/asw_dds.rds")

##create dds object for location_group analysis
dds_location_behaviour <- copy(dds)
##create groupings of weevil location and behavioural response to parasitoid
dds_location_behaviour$group <- factor(paste(dds$Weevil_Location,dds$Behaviour,sep="_"))

##add group to design
design(dds_location_behaviour) <- ~group
##run deseq2 and generate results
dds_location_behaviour <- DESeq(dds_location_behaviour)
##save dds_group
saveRDS(dds_location_behaviour, file = "output/deseq2/asw/WT_location_behaviour/asw_dds_location_behaviour.rds")