library("tximport")
library("data.table")
library("DESeq2")
library("VennDiagram")

sample_data <- fread("data/sample_key.csv")
trinotate <- fread("data/asw_transcriptome/most_sig_transcript_blastx_hit_for_each_gene.csv")

dds <- readRDS("output/deseq2/dds.rds")
##create dds object for group analysis
dds_group_parasitism <- copy(dds)

#######

##create groupings of weevil location and behavioural response to parasitoid
dds_group_parasitism$group <- factor(paste(dds$Weevil_Location,dds$Viral_expressed, sep="_"))

##add group to design
design(dds_group_parasitism) <- ~group
##run deseq2 and generate results
dds_group_parasitism <- DESeq(dds_group_parasitism)
##save dds_group
saveRDS(dds_group_parasitism, file = "output/deseq2/parasitism_status/dds_group_parasitism.rds")

resultsNames(dds_group_parasitism)

##Ru ev non-ev unpara: 0 DEGs
##Ru ev non-ev para: 5 DEGs (ribosomal, l2efl, no annots? - probably noise?)
##Dun ev non-ev para: 0 DEGs
##Dun ev non-ev unpara: 2 DEGs (transposon tigger, no annot - noise?)

##Ru viral vs non-viral: 1231 DEGs - would unpara ASW here express something resistancey that unpara Dun wouldn't?
##Dun viral vs non-viral: 1194 DEGs - compare (espl genes not de, nucleoprotein genes are)
##what's the overlap there?

###########
##Ruakura##
###########

rua_res_group <- results(dds_group_parasitism, contrast = c("group", "Ruakura_Yes", "Ruakura_No"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
rua_ordered_res_group <- rua_res_group[order(rua_res_group$padj),]
##Make data table and write to output
rua_ordered_res_group_table <- data.table(data.frame(rua_ordered_res_group), keep.rownames = TRUE)
rua_ordered_sig_res_group_table <- subset(rua_ordered_res_group_table, padj < 0.05)
rua_sig_annots <- merge(rua_ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE, all.y=FALSE)
##write full and sig res files
fwrite(rua_ordered_res_group_table, "output/deseq2/parasitism_status/ruakura_res_group.csv")
fwrite(rua_sig_annots, "output/deseq2/parasitism_status/ruakura_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)

###########
##Dunedin##
###########

dun_res_group <- results(dds_group_parasitism, contrast = c("group", "Dunedin_Yes", "Dunedin_No"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
dun_ordered_res_group <- dun_res_group[order(dun_res_group$padj),]
##Make data table and write to output
dun_ordered_res_group_table <- data.table(data.frame(dun_ordered_res_group), keep.rownames = TRUE)
dun_ordered_sig_res_group_table <- subset(dun_ordered_res_group_table, padj < 0.05)
dun_sig_annots <- merge(dun_ordered_sig_res_group_table, trinotate, by.x="rn", by.y="#gene_id", all.x=TRUE, all.y=FALSE)
##write full and sig res files
fwrite(dun_ordered_res_group_table, "output/deseq2/parasitism_status/dunedin_res_group.csv")
fwrite(dun_sig_annots, "output/deseq2/parasitism_status/dunedin_analysis_sig_degs.csv", col.names = TRUE, row.names = FALSE)

#######################
##Compare DEG Overlap##
#######################

#Draw Venn Diagram
vd <- venn.diagram(x = list("Ruakura"=rua_ordered_sig_res_group_table$rn, "Dunedin"=dun_ordered_sig_res_group_table$rn), filename=NULL, alpha=0.5, cex = 1, cat.cex=1, lwd=1)
grid.newpage()
grid.draw(vd)

##Ruakura specific
ru_specific_degs <- data.table(setdiff(rua_ordered_sig_res_group_table$rn, dun_ordered_sig_res_group_table$rn))
ru_specific_annots <- merge(rua_sig_annots, ru_specific_degs, by.x="rn", by.y="V1", all.x=FALSE, all.y=TRUE)

##Dun specific
dun_specific_degs <- data.table(setdiff(dun_ordered_sig_res_group_table$rn, rua_ordered_sig_res_group_table$rn))
dun_specific_annots <- merge(dun_sig_annots, dun_specific_degs, by.x="rn", by.y="V1", all.x=FALSE, all.y=TRUE)


##blastx for unann DEGs



