library(Biostrings)
library(tidyverse)

global_top3hit <- readDNAStringSet("../../data/misc/global_top3hit.fasta")
global_top3hit <- unique(global_top3hit)
# metadata_gisaid <- read_tsv("../../data/misc/metadata_2021-02-16_09-24.tsv")
formated_name_top3hit <- sapply(names(global_top3hit), function(x){
    gisaid_id_t <- strsplit(x, "|", fixed = T)[[1]][2]
    vn <- paste0("hCoV-19/", strsplit(x, "|", fixed = T)[[1]][1])
    date_t <- strsplit(x, "|", fixed = T)[[1]][3]
    location_t <- strsplit(vn, "/", fixed = T)[[1]][2]
    paste0(gisaid_id_t, "|", vn, "|", date_t, "|", location_t)
})
names(global_top3hit) <- formated_name_top3hit

gisaid_id_t <- sapply(names(global_top3hit), function(x){
    strsplit(x, "|", fixed = T)[[1]][2]
})
metadata_gisaid %>% filter(gisaid_epi_isl %in% gisaid_id_t) %>% write_csv("../../data/seq_global_analysis/metadata_global_top3hit.csv")

## deal with PANGO ref fasta
# EPI_ISL_416362|Shanghai/SH0054/2020|2020-02-02|Asia
seq_pango_ref <- readDNAStringSet("../../data/misc/unmasked_1275_pango_refs/1275_aligned_ref_29903.fasta")
seq_pango_ref <- seq_pango_ref[names(seq_pango_ref) != "MN908947_3"]

formated_name_pango_ref <- sapply(names(seq_pango_ref), function(x){
    gisaid_id_t <- strsplit(x, "|", fixed = T)[[1]][2]
    vn <- strsplit(x, "|", fixed = T)[[1]][1]
    date_t <- strsplit(x, "|", fixed = T)[[1]][3]
    location_t <- strsplit(vn, "/", fixed = T)[[1]][2]
    paste0(gisaid_id_t, "|", vn, "|", date_t, "|", location_t)
})
names(seq_pango_ref) <- formated_name_pango_ref

seq_global_all <- c(global_top3hit, seq_pango_ref)
names(seq_global_all) <- gsub(" ", "_", names(seq_global_all))

# mask gisaid masking position and head/tail 100
df_mask_gisaid <- read_csv("../../data/masked_sites_GISIAD.csv")
tmp <- rep(FALSE, width(seq_global_all)[1])
tmp[df_mask_gisaid$masked_sites] <- TRUE
tmp[1:100] <- TRUE
tmp[(width(seq_global_all)[1]-99):width(seq_global_all)[1]] <- TRUE

at <- matrix(rep(tmp, length(seq_global_all)),
            nrow=length(seq_global_all), ncol=width(seq_global_all)[1],
            byrow=TRUE)
letter_subject <- DNAString(paste(rep.int("-", width(seq_global_all)[1]), collapse=""))
letter <- as(Views(letter_subject, start=1, end=rowSums(at)), "XStringSet")
seq_global_all_GISIAD_masked <- replaceLetterAt(seq_global_all, at, letter)

gisaid_id_t <- sapply(names(seq_global_all_GISIAD_masked), function(x){
    strsplit(x, "|", fixed = T)[[1]][1]
})
seq_global_all_GISIAD_masked_rmdup <- seq_global_all_GISIAD_masked[!duplicated(gisaid_id_t)]

writeXStringSet(seq_global_all_GISIAD_masked, "../../data/seq_global_analysis//seq_global_all.fasta")
writeXStringSet(seq_global_all_GISIAD_masked_rmdup, "../../data/seq_global_analysis/seq_global_all_rmdup.fasta")

date_t <- sapply(names(seq_global_all_GISIAD_masked_rmdup), function(x){
    strsplit(x, "|", fixed = T)[[1]][3]
})
df_date <- tibble(seq_name = names(date_t), date = date_t)
write_tsv(df_date, "../../data/seq_global_analysis/date_global_all_rmdup.tsv", col_names = F)

