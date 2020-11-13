library("tximport")
library("data.table")
library("DESeq2")

dds <- readRDS("output/deseq2/mh/mh_dds.rds")
mh_trinotate_report <- fread("data/mh_edited_transcript_ids/trinotate_longest_isoform.csv")
##what about best hit per gene file?
mh_viral_blast_res <- fread("data/mh_edited_transcript_ids/mh_recip_viral_blast.csv")

##create dds object for location_group analysis
mh_dds_PCR_viral <- copy(dds)
##create groupings of weevil location and parasitism response to parasitoid
mh_dds_PCR_viral$group <- factor(paste(dds$Parasitism_status,dds$Viral_expressed,sep="_"))

##add group to design
design(mh_dds_PCR_viral) <- ~group
##run deseq2 and generate results
mh_dds_PCR_viral <- DESeq(mh_dds_PCR_viral)
resultsNames(mh_dds_PCR_viral)

##compare unpara-virus-expressing to unpara-not-virus-expressing
para_virus_res_group <- results(mh_dds_PCR_viral, contrast = c("group", "undetected_No", "undetected_Yes"), lfcThreshold = 1, alpha = 0.1)
##Order based of padj
para_virus_ordered_res_group <- para_virus_res_group[order(para_virus_res_group$padj),]
##Make data table and write to output
para_virus_ordered_res_group_table <- data.table(data.frame(para_virus_ordered_res_group), keep.rownames = TRUE)
para_virus_sig_res <- subset(para_virus_ordered_res_group_table, padj < 0.05)
para_virus_sig_res_annots <- merge(para_virus_sig_res, mh_trinotate_report, by.x="rn", by.y="#gene_id", all.x=TRUE)

mh_viral_blast_res$gene_id <- tstrsplit(mh_viral_blast_res$transcript_id, "_i", keep=c(1))
mh_viral_blast_res <- mh_viral_blast_res[,-c(1)]
mh_viral_blast_res$edited_gene_id <- paste("MH_",mh_viral_blast_res$gene_id, sep="")
para_virus_sig_recip_annots <- merge(para_virus_sig_res_annots, mh_viral_blast_res, by.x="rn", by.y="edited_gene_id", all.x=TRUE)
fwrite(para_virus_sig_recip_annots, "output/deseq2/mh/WT_unpara_viral/DEGs_annots.csv")

##they seem to DE many viral genes and RNAPs but not venom genes (unless hiding in those unann)

plotCounts(mh_dds_PCR_viral, "MH_TRINITY_DN13422_c0_g1", intgroup = c("group"), main="")
EnhancedVolcano(dun_viral_infection_sig_res, x="log2FoldChange", y="padj", lab="", pointSize = 3)
