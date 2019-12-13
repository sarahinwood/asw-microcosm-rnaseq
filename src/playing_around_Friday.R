dds <- readRDS("output/deseq2/dds.rds")
##create dds object for group analysis
dds_group <- copy(dds)

#######

##create groupings of weevil location and behavioural response to parasitoid
dds_group$Weevil_Location <- factor(dds$Weevil_Location)
dds_group$Behaviour <- factor(dds$Behaviour)
dds_group$Viral_expressed<- factor(dds$Viral_expressed)
##add group to design
design(dds_group) <- ~Weevil_Location+Viral_expressed+Behaviour
##run deseq2 and generate results
dds_group <- DESeq(dds_group)
##save dds_group
saveRDS(dds_group, file = "output/deseq2/dds_group.rds")


dds_abdo$Treatment <- factor(dds_abdo$Treatment, levels=time_order)
dds_abdo$Wasp_Location <- factor(dds_abdo$Wasp_Location)
dds_abdo$Flow_cell <- factor(dds_abdo$Flow_cell)
##add factors of ineterst to design
design(dds_abdo) <- ~Flow_cell+Wasp_Location+Treatment