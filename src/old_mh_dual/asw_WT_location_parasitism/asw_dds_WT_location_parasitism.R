library(data.table)
library(DESeq2)

##combine location + para status
##pairwise comp ru-undetected to dun-undetected
##(This doesn't consider expression in parasitised samples, or location-based differences)

##DEGs between unpara_ru and unpara_dun often high in all ru or dun (incl. para samples)
##not what we want when looking for resistance

dds <- readRDS("output/deseq2/asw/asw_dds.rds")

##create dds object for location_group analysis
dds_location_parasitism <- copy(dds)
##create groupings of weevil location and parasitism response to parasitoid
dds_location_parasitism$group <- factor(paste(dds$Weevil_Location,dds$Parasitism_status,sep="_"))

##add group to design
design(dds_location_parasitism) <- ~group
##run deseq2 and generate results
dds_location_parasitism <- DESeq(dds_location_parasitism)
##save dds_location_parasitism
saveRDS(dds_location_parasitism, file = "output/deseq2/asw/WT_location_parasitism/dds_location_parasitism.rds")

