library(data.table)
library(DESeq2)
library(ggplot2)
library(viridis)

##asw reads mapped
dds_concat_group_asw <- readRDS("output/deseq2/asw_dual/asw_dual_dds.rds")
counts_table_asw <- (data.table(counts(dds_concat_group_asw)))
counts_colSums_asw <- setDT(data.frame(colSums(counts_table_asw, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_asw, old=c("rn", "colSums.counts_table_asw..na.rm...TRUE."), new=c("sample_name", "readpairs_mapped_ASW"))

##mh reads mapped
dds_concat_group_mh <- readRDS("output/deseq2/mh_dual/mh_dual_dds.rds")
counts_table_mh <- (data.table(counts(dds_concat_group_mh)))
counts_colSums_mh <- setDT(data.frame(colSums(counts_table_mh, na.rm=TRUE)), keep.rownames=TRUE)
setnames(counts_colSums_mh, old=c("rn", "colSums.counts_table_mh..na.rm...TRUE."), new=c("sample_name", "readpairs_mapped_MH"))

##both reads mapped
read_mapping <- merge(counts_colSums_asw, counts_colSums_mh)
bbduk_reads_out <- fread("output/bbduk_trim/bbduk_reads_out.csv")
full_read_mapping <- merge(read_mapping, bbduk_reads_out, by="sample_name")
salmon_mapping <- fread("output/asw_mh_concat_salmon/salmon_mapping.csv")
full_read_mapping <- merge(full_read_mapping, salmon_mapping)

##mapping %s
full_read_mapping$total_mapped_readpairs <- (full_read_mapping$readpairs_mapped_MH + full_read_mapping$readpairs_mapped_ASW)
full_read_mapping$`%_ofmapped_ASW` <- (full_read_mapping$readpairs_mapped_ASW/full_read_mapping$total_mapped_readpairs)*100
full_read_mapping$`%_ofmapped_MH` <- (full_read_mapping$readpairs_mapped_MH/full_read_mapping$total_mapped_readpairs)*100
fwrite(full_read_mapping, "output/deseq2/read_mapping.csv")

para <- subset(full_read_mapping, Parasitism_status=="parasitized")
mean(para$`%_ofmapped_MH`)
undetected <- subset(full_read_mapping, Parasitism_status=="undetected")
mean(undetected$`%_ofmapped_MH`)

#full_read_mapping$para_loc <- paste(full_read_mapping$Parasitism_status, full_read_mapping$Location, sep="_")

#dun<-subset(full_read_mapping, full_read_mapping$Location == "Dunedin")

##plot or calc means for mh mapping between para and unpara
ggplot(full_read_mapping, aes(x=Parasitism_status, y=`%_ofmapped_MH`))+
  geom_boxplot(aes(fill=Parasitism_status, colour=Parasitism_status), alpha=0.7, show.legend = FALSE)+
  theme_bw()+
  scale_fill_viridis(discrete=TRUE)+
  scale_colour_viridis(discrete=TRUE)+
  xlab("Parasitism status")+
  ylab(expression(paste("% of normalized counts to ", italic("M. hyperodae"), " transcriptome")))
  
##t test for Mh mapping between para and undetected
t.test(`%_ofmapped_MH` ~ Parasitism_status, data=full_read_mapping)

#t test for mapping rates of parasitised ASW on location
para_mapping <- subset(full_read_mapping, Parasitism_status=="parasitized")
para_mapping$location <- tstrsplit(para_mapping$sample_name, "", keep=c(1))
t.test(`%_ofmapped_MH` ~location, data=para_mapping)

##plot or calc means for mh mapping between para and unpara
ggplot(full_read_mapping, aes(x=Parasitism_status, y=`%_ofmapped_ASW`))+
  geom_boxplot(aes(fill=Parasitism_status, colour=Parasitism_status), alpha=0.7, show.legend = FALSE)+
  theme_bw()+
  scale_fill_viridis(discrete=TRUE)+
  scale_colour_viridis(discrete=TRUE)+
  xlab("Parasitism status")+
  ylab(expression(paste("% of normalized counts to ASW transcriptome")))
