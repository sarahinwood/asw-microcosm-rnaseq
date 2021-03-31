library("data.table")
library("DESeq2")

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
##save dds_location_behaviour
saveRDS(dds_location_behaviour, file = "output/deseq2/asw/WT_location_behaviour/asw_dds_location_behaviour.rds")
dds_location_behaviour <- readRDS("output/deseq2/asw/WT_location_behaviour/asw_dds_location_behaviour.rds")

###########
##Ruakura##
###########

rua_res_group <- results(dds_location_behaviour, contrast = c("group", "Ruakura_E", "Ruakura_N"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
rua_ordered_res_group <- rua_res_group[order(rua_res_group$padj),]
##Make data table and write to output
rua_ordered_res_group_table <- data.table(data.frame(rua_ordered_res_group), keep.rownames = TRUE)
rua_ordered_sig_res_group_table <- subset(rua_ordered_res_group_table, padj < 0.05)

##write full and sig res files
fwrite(rua_ordered_res_group_table, "output/deseq2/asw/WT_location_behaviour/ruakura_res_group.csv")
fwrite(rua_ordered_sig_res_group_table, "output/deseq2/asw/WT_location_behaviour/ruakura_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)

###########
##Dunedin##
###########

dun_res_group <- results(dds_location_behaviour, contrast = c("group", "Dunedin_E", "Dunedin_N"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
dun_ordered_res_group <- dun_res_group[order(dun_res_group$padj),]
##Make data table and write to output
dun_ordered_res_group_table <- data.table(data.frame(dun_ordered_res_group), keep.rownames = TRUE)
dun_ordered_sig_res_group_table <- subset(dun_ordered_res_group_table, padj < 0.05)

##write full and sig res files
fwrite(dun_ordered_res_group_table, "output/deseq2/asw/WT_location_behaviour/dunedin_res_group.csv")
fwrite(dun_ordered_sig_res_group_table, "output/deseq2/asw/WT_location_behaviour/dunedin_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)


