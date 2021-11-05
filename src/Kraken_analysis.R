library(data.table)
library(tidyverse)

############################
## Pavian report analysis ##
############################

pavian_res <- fread("output/kraken/pavian/pavian_results_overview.csv")
pavian_res$sample_name <- tstrsplit(pavian_res$Name, "_", keep=c(2))

sample_table <- fread("data/sample_table.csv")
sample_info <- sample_table[,c(1,9,12,14,15)]

pavian_sample_info <- merge(pavian_res, sample_info)
##remove sample with 12% classified
#pavian_sample_info <- subset(pavian_sample_info, pavian_sample_info$sample_name!="RY14C")

##test normality of bacterial reads
shapiro.test(pavian_sample_info$pct_bacterial)
shapiro.test(pavian_sample_info$pct_viral)
##distribution non-normal - but have 30+ samples so think its okay to use parametric still

## BACTERIAL ##
##bacterial location - not sig
t.test(pct_bacterial ~ Weevil_Location, data = pavian_sample_info)
##bacterial attack - not sig
t.test(pct_bacterial ~ Attacked, data = pavian_sample_info)
##bacterial para
        ##sig diff in bacterial between para and undetected
        ## higher bacterial in unpara (3.1% vs 2.5% in para)
t.test(pct_bacterial ~ Parasitism_status, data = pavian_sample_info)

## VIRAL ##
##viral location - not sig
t.test(pct_viral ~ Weevil_Location, data = pavian_sample_info)
##viral attack - not sig
t.test(pct_viral ~ Attacked, data = pavian_sample_info)
##viral para - sig diff in bacterial between para and undetected? - not sig
t.test(pct_viral ~ Parasitism_status, data = pavian_sample_info)
##viral vs viral expression - not sig - likely doesn't classify MhV reads as LbFV pretty weird and distinct
t.test(pct_viral ~ Viral_expressed, data = pavian_sample_info)

#########################
## full kraken results ##
#########################
##what are bacterial reads in para and unpara given sig diff?

# read file path
files <- list.files(path = "output/kraken/reports", pattern = "*.txt", full.names = TRUE)
##read in
kraken_res <- rbindlist(lapply(files, fread), idcol = "origin")
##colname for file name
kraken_res[, origin := factor(origin, labels = basename(files))]
##tidy fileame to sample name
kraken_res$sample_name <- tstrsplit(kraken_res$origin, "_", keep=c(2))

setnames(kraken_res, old=c("V1", "V2", "V3", "V4", "V5", "V6"),
         new=c("percent_classified", "number_clade", "number_taxon", "rank", "taxid", "scientific_name"))
kraken_res_table <- kraken_res[,c(8, 2, 5,6,7)]
##remove rows with 0%
#kraken_res_table_small <- subset(kraken_res_table, percent_classified>0)

##down to phyla level
kraken_res_table_phylum <- subset(kraken_res_table, rank=="P")
phlya_summary <- aggregate(percent_classified ~ scientific_name, data=kraken_res_table_phylum, mean)

##taxid list of bacterial taxa?
taxid_bacterial <- fread("output/taxids/bacterial_taxids.txt", header=FALSE)
kraken_res_table_phylum_bacterial <- subset(kraken_res_table_phylum, taxid %in% taxid_bacterial$V1)
##then keep bacterial classes where a sample has above 0 counts
kraken_res_table_phylum_bacterial <- subset(kraken_res_table_phylum_bacterial, percent_classified>0)
bacteria_phlya_summary <- aggregate(percent_classified ~ scientific_name, data=kraken_res_table_phylum_bacterial, mean)

##most reads are proteobacteria
proteobacteria <- subset(kraken_res_table, kraken_res_table$scientific_name=="Proteobacteria")
proteobacteria_info <- merge(proteobacteria, sample_info)
##para
t.test(percent_classified ~ Parasitism_status, data = proteobacteria_info)
##attack
t.test(percent_classified ~ Attacked, data = proteobacteria_info)
##location
t.test(percent_classified ~ Weevil_Location, data = proteobacteria_info)

##class level
kraken_res_table_class <- subset(kraken_res_table, rank=="C")
kraken_res_table_class_bacterial <- subset(kraken_res_table_class, taxid %in% taxid_bacterial$V1)
bacterial_class_summary <- aggregate(percent_classified ~ scientific_name, data=kraken_res_table_class_bacterial, mean)

##order level
kraken_res_table_order <- subset(kraken_res_table, rank=="O")
kraken_res_table_order_bacterial <- subset(kraken_res_table_order, taxid %in% taxid_bacterial$V1)
bacterial_order_summary <- aggregate(percent_classified ~ scientific_name, data=kraken_res_table_order_bacterial, mean)

Enterobacterales <- subset(kraken_res_table_order_bacterial, scientific_name=="Enterobacterales")
Enterobacterales_info <- merge(Enterobacterales, sample_info)
t.test(percent_classified ~ Parasitism_status, data = Enterobacterales_info)

Xanthomonadales <- subset(kraken_res_table_order_bacterial, scientific_name=="Xanthomonadales")
Xanthomonadales_info <- merge(Xanthomonadales, sample_info)
t.test(percent_classified ~ Parasitism_status, data = Xanthomonadales_info)

Flavobacteriales <- subset(kraken_res_table_order_bacterial, scientific_name=="Flavobacteriales")
Flavobacteriales_info <- merge(Flavobacteriales, sample_info)
t.test(percent_classified ~ Parasitism_status, data = Flavobacteriales_info)
