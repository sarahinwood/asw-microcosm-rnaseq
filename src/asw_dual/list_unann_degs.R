library(data.table)
library(dplyr)

##read in all lists of DEGs
para <- fread("output/deseq2/asw_dual/WT_parasitism/sig_annots.csv")
location <- fread("output/deseq2/asw_dual/WT_location/sig_annots.csv")
attacked <- fread("output/deseq2/asw_dual/WT_attacked/sig_annots.csv")
para_loc <- fread("output/deseq2/asw_dual/INT_WT_parasitism-location/sig_annots.csv")
att_loc <- fread("output/deseq2/asw_dual/INT_WT_attacked-location/sig_annots.csv")

##make full table of all DEGs
full_table <- full_join(para, location)
full_table <- full_join(full_table, attacked)
full_table <- full_join(full_table, para_loc)
full_table <- full_join(full_table, att_loc)
##write list for blast
sig_no_annots <- subset(full_table, full_table$sprot_Top_BLASTX_hit=="")
fwrite(list(sig_no_annots$rn), "output/deseq2/asw_dual/no_annot/DEGs_ID_no_annot.txt")

##blastx for unann
unann_deg_blast <- fread("output/deseq2/asw_dual/no_annot/blastx.outfmt6")
setnames(unann_deg_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("transcript_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
unann_annots <- unann_deg_blast[,c(1,11,13)]
unann_annots$rn <- tstrsplit(unann_annots$transcript_id, "_i", keep=c(1))
##merge unann blast with other DEGs
##parasitism
para_annots <- merge(para, unann_annots, by="rn", all.x=TRUE)
fwrite(para_annots, "output/deseq2/asw_dual/WT_parasitism/sig_blast_annots.csv")
##location
location_annots <- merge(location, unann_annots, by="rn", all.x=TRUE)
fwrite(location_annots, "output/deseq2/asw_dual/WT_location/sig_blast_annots.csv")
##attack
attacked_annots <- merge(attacked, unann_annots, by="rn", all.x=TRUE)
fwrite(attacked_annots, "output/deseq2/asw_dual/WT_attacked/sig_blast_annots.csv")
##para-loc
para_loc_annots <- merge(para_loc, unann_annots, by="rn", all.x=TRUE)
fwrite(para_loc_annots, "output/deseq2/asw_dual/INT_WT_parasitism-location/sig_blast_annots.csv")
##attack-loc
att_loc_annots <- merge(att_loc, unann_annots, by="rn", all.x=TRUE)
fwrite(att_loc_annots, "output/deseq2/asw_dual/INT_WT_attacked-location/sig_blast_annots.csv")
