library(data.table)
library(DESeq2)

##combine location + para status
##pairwise comp ru-undetected to dun-undetected
##(This doesn't consider expression in parasitised samples, or location-based differences)

##DEGs between unpara_ru and unpara_dun often high in all ru or dun (incl. para samples)
##not what we want when looking for resistance

dds <- readRDS("output/deseq2/asw/asw_dds.rds")

##create dds object for location_group analysis
dds_location_attacked <- copy(dds)
##create groupings of weevil location and attacked response to parasitoid
dds_location_attacked$group <- factor(paste(dds$Weevil_Location,dds$Attacked,sep="_"))

##add group to design
design(dds_location_attacked) <- ~group
##run deseq2 and generate results
dds_location_attacked <- DESeq(dds_location_attacked)
##save dds_location_attacked
saveRDS(dds_location_attacked, file = "output/deseq2/asw/WT_location_attacked/dds_location_attacked.rds")

