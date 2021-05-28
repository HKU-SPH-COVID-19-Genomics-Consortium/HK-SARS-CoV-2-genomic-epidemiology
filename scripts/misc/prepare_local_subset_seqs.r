library(Biostrings)
library(tidyverse)

seqs_local <- readDNAStringSet("../../data/seq_global_analysis/hk_case_primer_masked_GISAID_masked_2021-05-08.afa")
seqs_local <- seqs_local[names(seqs_local) != "MN908947_3"]
case_id <- sapply(names(seqs_local), function(x){
	tmp <- strsplit(x, "|", fixed = T)[[1]]
	tmp <- tmp[grepl("HKcase_", tmp)]
	strsplit(tmp, "_")[[1]][2]
})

w3_all <- read_tsv("../../data/misc/w3_location.txt")
w3_need <- read_tsv("../../data/misc/w3_sub_100_ids.txt", col_names = F)
w4_all <- read_tsv("../../data/misc/w4_location.txt")
w4_need <- read_tsv("../../data/misc/w4_sub_65_ids.txt", col_names = F)

w3_rm <- w3_all$traits[!w3_all$traits %in% w3_need$X1]
w4_rm <- w4_all$traits[!w4_all$traits %in% w4_need$X1]

to_rm <- c(w3_rm, w4_rm)
to_rm <- sapply(to_rm, function(x) {
	tmp <- strsplit(x, "/", fixed = T)[[1]][1]
	tmp <- strsplit(tmp, "_", fixed = T)[[1]][2]
})

seqs_local_subset <- seqs_local[!case_id %in% to_rm]
writeXStringSet(seqs_local_subset, "../../data/seq_global_analysis/hk_case_primer_masked_GISAID_masked_2021-05-08_subset.afa")
