library(data.table)
library(dplyr)

##read in all lists of DEGs
para <- fread("output/deseq2/asw_dual/WT_parasitism/sig_annots.csv")
location <- fread("output/deseq2/asw_dual/WT_location/sig_annots.csv")
attacked <- fread("output/deseq2/asw_dual/WT_attacked/sig_annots.csv")
para_loc <- fread("output/deseq2/asw_dual/INT_WT_parasitism-location/sig_annots.csv")
att_loc <- fread("output/deseq2/asw_dual/INT_WT_attacked-location/sig_annots.csv")

##make full table of all DEGs
full_table <- full_join(para, full_join(location, full_join(attacked, full_join(para_loc, att_loc))))

##write list for blast
sig_no_annots <- subset(full_table, full_table$sprot_Top_BLASTX_hit=="")
sig_no_annots$right_rn <- tstrsplit(sig_no_annots$rn, ".", fixed=TRUE, keep=c(1)) 
fwrite(list(unique(sig_no_annots$right_rn)), "output/deseq2/asw_dual/no_annot/DEGs_ID_no_annot.txt")

##blastx for unann
unann_deg_blast <- fread("output/deseq2/asw_dual/no_annot/blastx.outfmt6")
setnames(unann_deg_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("transcript_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##remove unchar/hypo annots
no_hypo <- subset(unann_deg_blast, !grepl("hypothetical|unchar|GSCOCG|unnamed protein product", annotation))
##sort for best hit
setorder(no_hypo, transcript_id, evalue, -bit_score)
min_evalues <- no_hypo[,.SD[which.min(evalue)], by=transcript_id]

##best hit from all
setorder(unann_deg_blast, transcript_id, evalue, -bit_score)
min_hypo <- unann_deg_blast[,.SD[which.min(evalue)], by=transcript_id]
only_hypo <- setdiff(min_hypo$transcript_id, min_evalues$transcript_id)
min_hypo_only <- subset(min_hypo, transcript_id %in% only_hypo)
##ASW_TRINITY_DN48834_c0_g1_i1 - pieris annot - like MhV? best hit to virus - attacked
##compare to MhV pieris hits?

##ASW_TRINITY_DN8033_c0_g1_i1 best hit wuhan virus? - para loc

min_evalues$rn <- tstrsplit(min_evalues$transcript_id, "_i", keep=c(1))
##merge unann blast with other DEGs
##parasitism
para_annots <- merge(para, min_evalues, by="rn", all.x=TRUE)
fwrite(para_annots, "output/deseq2/asw_dual/WT_parasitism/sig_blast_annots.csv")
##location
location_annots <- merge(location, min_evalues, by="rn", all.x=TRUE)
fwrite(location_annots, "output/deseq2/asw_dual/WT_location/sig_blast_annots.csv")
##attack
attacked_annots <- merge(attacked, min_evalues, by="rn", all.x=TRUE)
fwrite(attacked_annots, "output/deseq2/asw_dual/WT_attacked/sig_blast_annots.csv")
##para-loc
para_loc_annots <- merge(para_loc, min_evalues, by="rn", all.x=TRUE)
fwrite(para_loc_annots, "output/deseq2/asw_dual/INT_WT_parasitism-location/sig_blast_annots.csv")
##attack-loc
att_loc_annots <- merge(att_loc, min_evalues, by="rn", all.x=TRUE)
fwrite(att_loc_annots, "output/deseq2/asw_dual/INT_WT_attacked-location/sig_blast_annots.csv")
