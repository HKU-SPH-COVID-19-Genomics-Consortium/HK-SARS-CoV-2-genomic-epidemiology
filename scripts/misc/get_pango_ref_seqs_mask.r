library(tidyverse)
library(Biostrings)

metadata_gisaid <- read_delim("../../data/metadata_2021-02-16_09-24.tsv", "\t", , col_types = cols(.default = "c"))
gisaid_fasta <- readDNAStringSet("../../data/mmsa_2021-02-16.fa")
# gisaid_fasta_m <- consensusMatrix(gisaid_fasta)
list_pango_ref <- readLines("../../data/nextstrain_pango_ref.txt")
id_pango_ref <- sapply(list_pango_ref, function(x){
	strsplit(x, "|", fixed = T)[[1]][1]
})
seqs_pango_ref <- gisaid_fasta[names(gisaid_fasta) %in% id_pango_ref]

formated_name <- sapply(names(seqs_pango_ref), function(x){
    df_i <- metadata_gisaid %>% filter(gisaid_epi_isl == x)
    paste0(df_i$gisaid_epi_isl, "|", df_i$strain, "|", df_i$date, "|", df_i$region)
})
names(seqs_pango_ref) <- formated_name

# mask gisaid masking position
df_mask_gisaid <- read_csv("../../data/masked_sites_GISIAD.csv")
tmp <- rep(FALSE, width(seqs_pango_ref)[1])
tmp[df_mask_gisaid$masked_sites] <- TRUE
sum(!df_mask_gisaid$masked_sites %in% c(1:100, (29903-99):29903))

at <- matrix(rep(tmp, length(seqs_pango_ref)),
            nrow=length(seqs_pango_ref), ncol=width(seqs_pango_ref)[1],
            byrow=TRUE)
letter_subject <- DNAString(paste(rep.int("-", width(seqs_pango_ref)[1]), collapse=""))
letter <- as(Views(letter_subject, start=1, end=rowSums(at)), "XStringSet")
seqs_pango_ref_GISIAD_masked <- replaceLetterAt(seqs_pango_ref, at, letter)
writeXStringSet(seqs_pango_ref_GISIAD_masked, "../../results/pango_ref.fasta")