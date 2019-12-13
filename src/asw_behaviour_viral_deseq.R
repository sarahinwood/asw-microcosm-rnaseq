library(tximport)
library(data.table)
library(DESeq2)
library(dplyr)

##Import table describing samples
sample_data <- fread("data/sample_key.csv", header=TRUE)
setkey(sample_data, Sample_name)

##read in dds from non-viral script
dds <- readRDS("output/deseq2/dds.rds")
##create dds object for group analysis
dds_viral <- copy(dds)

##create groupings of weevil location and behavioural response to parasitoid
dds_viral$group <- factor(paste(dds$Weevil_Location,dds$Behaviour,dds$Viral_expressed,sep="_"))

##add group to design
design(dds_viral) <- ~group
##run deseq2 and generate results
dds_viral <- DESeq(dds_viral)
##save dds_viral
saveRDS(dds_viral, file = "output/deseq2/dds_viral.rds")

##look at raw and normalised counts for viral genes

plotCounts(dds_viral, "TRINITY_DN38122_c0_g1", intgroup = c("group"), main = "Bro1")

#bro1 - TRINITY_DN38122_c0_g1
#bro2 - TRINITY_DN35519_c0_g1
#KilA-N - TRINITY_DN1684_c0_g1

##raw counts
counts_matrix_dt <- data.table(counts(dds_viral), keep.rownames = TRUE)
viral_genes <- list("TRINITY_DN38122_c0_g1", "TRINITY_DN35519_c0_g1", "TRINITY_DN1684_c0_g1")
viral_counts_matrix <- counts_matrix_dt[counts_matrix_dt$rn %in% viral_genes]
viral_counts_melted <- melt(viral_counts_matrix)

##get dt with normalised counts for each of the 3 viral genes and prep to merge with raw count table
kilA_N_normalized_counts <- setDT(plotCounts(dds_viral, "TRINITY_DN1684_c0_g1", intgroup = c("group"), returnData=TRUE), keep.rownames=TRUE)
kilA_N_normalized_counts$gene <- paste("norm_TRINITY_DN1684_c0_g1")
kilA_N_normalized_counts <- kilA_N_normalized_counts[, !"group"]
kilA_N_wide <- dcast(kilA_N_normalized_counts, gene ~ rn, value.var="count")

bro1_normalized_counts <- setDT(plotCounts(dds_viral, "TRINITY_DN38122_c0_g1", intgroup = c("group"), returnData=TRUE), keep.rownames=TRUE)
bro1_normalized_counts$gene <- paste("norm_TRINITY_DN38122_c0_g1")
bro1_normalized_counts <- bro1_normalized_counts[, !"group"]
bro1_wide <- dcast(bro1_normalized_counts, gene ~ rn, value.var="count")

bro2_normalized_counts <- setDT(plotCounts(dds_viral, "TRINITY_DN35519_c0_g1", intgroup = c("group"), returnData=TRUE), keep.rownames=TRUE)
bro2_normalized_counts$gene <- paste("norm_TRINITY_DN35519_c0_g1")
bro2_normalized_counts <- bro2_normalized_counts[, !"group"]
bro2_wide <- dcast(bro2_normalized_counts, gene ~ rn, value.var="count")

normalized_counts <- full_join(kilA_N_normalized_counts, bro1_normalized_counts)
normalized_counts <- data.table(full_join(normalized_counts, bro2_normalized_counts))
normalized_counts <- dcast(normalized_counts, gene~rn, value.var="count")
##merge raw and normalized counts
norm_raw_counts <- full_join(normalized_counts, viral_counts_matrix)
fwrite(norm_raw_counts, "output/deseq2/raw+normalized_counts.csv")
