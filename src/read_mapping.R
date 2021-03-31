library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)

##asw reads mapped
dds_concat_group_asw <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
counts_table_asw <- (data.table(counts(dds_concat_group_asw)))
counts_colSums_asw <- setDT(data.frame(colSums(counts_table_asw, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_asw, old=c("rn", "colSums.counts_table_asw..na.rm...TRUE."), new=c("Sample_name", "readpairs_mapped_ASW"))

##mh reads mapped
dds_concat_group_mh <- readRDS("output/deseq2/mh_dual/mh_dual_dds.rds")
counts_table_mh <- (data.table(counts(dds_concat_group_mh)))
counts_colSums_mh <- setDT(data.frame(colSums(counts_table_mh, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_mh, old=c("rn", "colSums.counts_table_mh..na.rm...TRUE."), new=c("Sample_name", "readpairs_mapped_MH"))

##both reads mapped
read_mapping <- merge(counts_colSums_asw, counts_colSums_mh)
bbduk_reads_out <- fread("output/bbduk_trim/bbduk_reads_out.csv")
full_read_mapping <- merge(read_mapping, bbduk_reads_out, by="Sample_name")
full_read_mapping$bbduk_halved <- (full_read_mapping$bbduk_reads_out)/2

##mapping %s
full_read_mapping$total_mapped_reads <- (full_read_mapping$readpairs_mapped_MH + full_read_mapping$readpairs_mapped_ASW)
full_read_mapping$`overall_mapping_%` <- (full_read_mapping$total_mapped_reads)/(full_read_mapping$bbduk_halved)*100
full_read_mapping$`%_ofmapped_ASW` <- (full_read_mapping$readpairs_mapped_ASW/full_read_mapping$total_mapped_reads)*100
full_read_mapping$`%_ofmapped_MH` <- (full_read_mapping$readpairs_mapped_MH/full_read_mapping$total_mapped_reads)*100
full_read_mapping$`%_ofall_ASW` <- (full_read_mapping$readpairs_mapped_ASW/full_read_mapping$bbduk_halved)*100
full_read_mapping$`%_ofall_MH` <- (full_read_mapping$readpairs_mapped_MH/full_read_mapping$bbduk_halved)*100

fwrite(full_read_mapping, "output/deseq2/read_mapping.csv")

full_read_mapping$para_loc <- paste(full_read_mapping$Parasitism_status, full_read_mapping$Location, sep="_")

dun<-subset(full_read_mapping, full_read_mapping$Location == "Dunedin")

##plot or calc means for mh mapping between para and unpara
ggplot(full_read_mapping, aes(x=Parasitism_status, y=`%_ofall_MH`))+
  geom_boxplot(aes(fill=Parasitism_status), outlier.shape=NA, alpha=0.8, show.legend = FALSE)+
  theme_bw()+
  scale_fill_viridis_d()+
  coord_cartesian(ylim=c(0, 10))+
  xlab("Parasitism status")+
  ylab(expression(paste("% reads mapping to ", italic("M. hyperodae"), " transcriptome")))
  
