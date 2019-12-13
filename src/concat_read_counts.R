library(data.table)
library(DESeq2)

dds_concat_group_asw <- readRDS("output/concat_deseq2/dds_concat_group_asw.rds")
counts_table_asw <- (data.table(counts(dds_concat_group_asw)))
counts_colSums_asw <- setDT(data.frame(colSums(counts_table_asw, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_asw, old=c("rn", "colSums.counts_table_asw..na.rm...TRUE."), new=c("Sample_name", "reads_mapped_ASW"))
##merge with csv with total reads
total_reads <- fread("output/concat_deseq2/sample_total_reads.csv")
total_reads_and_counts <- merge(counts_colSums_asw, total_reads,by="Sample_name")

dds_concat_group_mh <- readRDS("output/concat_deseq2/dds_concat_group_mh.rds")
counts_table_mh <- (data.table(counts(dds_concat_group_mh)))
counts_colSums_mh <- setDT(data.frame(colSums(counts_table_mh, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_mh, old=c("rn", "colSums.counts_table_mh..na.rm...TRUE."), new=c("Sample_name", "reads_mapped_MH"))
##merge with existing total reads and counts table
total_reads <- fread("output/concat_deseq2/sample_total_reads.csv")
total_reads_and_counts <- merge(counts_colSums_mh, total_reads_and_counts,by="Sample_name")
##calc. total mapped reads
total_reads_and_counts$total_mapped_reads <- (total_reads_and_counts$reads_mapped_MH + total_reads_and_counts$reads_mapped_ASW)
total_reads_and_counts$`%_ASW` <- (total_reads_and_counts$reads_mapped_ASW/total_reads_and_counts$total_mapped_reads)
total_reads_and_counts$`%_MH` <- (total_reads_and_counts$reads_mapped_MH/total_reads_and_counts$total_mapped_reads)
fwrite(total_reads_and_counts, "output/concat_deseq2/total_reads_and_counts.csv")


