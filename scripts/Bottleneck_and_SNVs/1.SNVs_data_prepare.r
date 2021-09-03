library(ggsci)
library(tidyverse)
library(ape)
library(parallel)
library(Biostrings)
library(seqinr)
library('RColorBrewer')
library(parallel)
library(ggrepel)
library(ggsignif)

# mc.cores=12
mc.cores=48

hk_seqs <- readDNAStringSet("../../data/seq_local_analysis/hk_case_primer_masked_2021-05-08.afa")
hk_seqs_1623 <- readDNAStringSet("../../data/hk_seq_1623.fasta")
seqed_meta <- read_csv("../../data/cov_and_lin.csv") 
id_case_df <- read_csv("../../data/labID_caseID_merged.csv", col_types = cols(.default = "c"))
data_metadata <- readxl::read_excel("../../data/Copy of 2019ncov-global-surveillance-linelist (2021.01.26).xlsx", col_types = "text")
data_metadata <- data_metadata %>% select(`HK case no.`:`onset to report (days)`)
data_metadata <- data_metadata[!is.na(data_metadata$`HK case no.`),]
data_metadata$case_id <- seq_len(nrow(data_metadata))
data_metadata$club_sim <- gsub("[^[:alnum:]]", "_", data_metadata$`Bar/pub/club`)
data_metadata <- data_metadata[!duplicated(data_metadata$`HK case no.`),]
metadata <- data_metadata
id_case_df_hk <- id_case_df %>% filter(case_ID %in% names(hk_seqs))
data_hk_study_samples <- read_csv("../../results/samples_hk_study.csv")
df_sample_type <- read_csv("../../results/wave_type_rp.csv")
data_duplicated_samples <- read_csv("../../results/duplicated_cases.csv")
data_duplicated_samples_keep <- data_duplicated_samples %>% filter(!(same_sample & diff_in_consensus>0)) 
data_duplicated_samples_keep <- data_duplicated_samples_keep %>% filter(`Genome coverage`*29903>=27000)
data_duplicated_samples_keep <- data_duplicated_samples_keep[data_duplicated_samples_keep$case_id %in% data_duplicated_samples_keep$case_id[duplicated(data_duplicated_samples_keep$case_id)],]
data_duplicated_samples_keep <- data_duplicated_samples_keep %>% filter(!case_id == "7409")
df_pair <- read_csv("../../results/Bottleneck_and_SNVs/bottleneck_1day.csv")
annot_gene <- read_tsv("../../data/SARS-CoV-2.gtf", col_names = F)

# clean annotation
annot_gene$Gene <- sapply(annot_gene$X9, function(x){
    strsplit(x, '"', fixed = T)[[1]][2]
})
unique(annot_gene$Gene)
annot_gene <- annot_gene %>% filter(Gene %in% c(paste0("nsp", 1:10), paste0("nsp", 12:16), "nsp12_1", "nsp12_2", "S", "ORF3a", "E", "M", "ORF6", "ORF7a_NOL", "ORF7b_NOL", "ORF8", "N", "ORF10")) %>% arrange(X4) %>% mutate(X9 = Gene) %>% select(-Gene)
annot_gene$X9 <- gsub("_NOL", "", annot_gene$X9) 
write_tsv(annot_gene, "../results/gene_annot.gtf")
annot_gene_sim <- annot_gene %>% select(X9, X4, X5)
names(annot_gene_sim) <- c("sequence","start","stop")
write_csv(annot_gene_sim, "../data/ORF_SCoV2.csv")

# read NGS readcount data
files <- list.files("../../2020-09-01_COVID_NGS_pipeline/results/", "data_snvs", full.names = T)
sample_files <- sapply(files, function(x){
    num <- strsplit(x, "/")[[1]]
    num <- num[length(num)]
    num <- strsplit(num, "-")[[1]][1]
    num <- strsplit(num, "_")[[1]]
    num <- num[length(num)]
    sample_tmp <- strsplit(num, ".csv")[[1]][1]
})
sum(sample_files %in% seqed_meta$Sample)

data_ori <- lapply(seq_along(files), function(i){
    x <- files[i]
    sample_tmp <- sample_files[i]
    case_id_tmp <- seqed_meta$`HK case no.`[seqed_meta$Sample == sample_tmp][1]
    print(case_id_tmp)
    if(length(case_id_tmp)==0){
        return(NA)
    }
    if(is.na(case_id_tmp)){
        return(NA)
    }
    check <- (sample_tmp %in% data_hk_study_samples$sample) | (sample_tmp %in% data_duplicated_samples_keep$`Lab No.`)
    # check <- TRUE
    if(check){ ## only include cases in the study
        tmp <- read_csv(x, col_types = cols(.default = "c"))
        # sample_t <- tmp$Sample[1]
        # sample_t <- strsplit(sample_t, "-")[[1]][1]
        # case_id_t <- id_case_df$case_ID[id_case_df$lab_ID == sample_t | id_case_df$lab_ID == tmp$Sample[1]][1]
        tmp$case_id <- case_id_tmp
        check0 <- length(unique(unlist(strsplit(tmp$`Identified by`, ", "))))>2
        # check1 <- seqed_meta$`Genome coverage`[seqed_meta$Sample == sample_tmp]*29903>=27000
        # check2 <- case_id_tmp %in% names(hk_seqs)
        if(check0){
            return(tmp)
        } else {
            return(NA)
        }
    } else {
        return(NA)
    }
})
data_ori <- data_ori[!is.na(data_ori)]
data_ori <- bind_rows(data_ori)
data_ori$case_id <- as.numeric(data_ori$case_id)
data_ori <- data_ori[!is.na(data_ori$case_id),]

unique(data_ori$Gene)

data_ori$Gene[is.na(data_ori$Gene)] <- "Not in ORF"
data_ori$Gene <- factor(data_ori$Gene, unique(data_ori$Gene))
data_ori$Silent_mutation[is.na(data_ori$Silent_mutation)] <- "Not in ORF"
data_ori <- data_ori %>% filter(!grepl("-S", Sample))
length(unique(data_ori$Sample)) 
length(unique(data_ori$case_id))

sample_should_have <- c(data_hk_study_samples$sample, data_duplicated_samples_keep$`Lab No.`)

# filter
# for the SNPs study, we only consider SNPs with depth >= 100, and was identified by more than one variant caller.
data_ori <- data_ori %>% filter(as.numeric(depth)>=100) %>% filter(grepl(",", `Identified by`, fixed = T))
data_ori$depth <- as.numeric(data_ori$depth)
data_ori$Alle_freq <- as.numeric(data_ori$Alle_freq)
data_ori$ID <- paste0(data_ori$case_id, "_", data_ori$Sample)
data_ori$mutation <- paste0(data_ori$`Ref base`, data_ori$Position, data_ori$`Alt base`)

# primer related masked sites
masked_sites <- read_csv("../../data/Bottleneck_and_SNVs/primer_masked_sites.csv")

# masking sites
masked_sites <- c(data_high_freq$Position, 1:100, (29903-99):29903)
data_ori <- data_ori %>% filter(!Position %in% masked_sites)

# label iSNVs
total_i <- nrow(data_ori)
data_ori$iSNVs <- sapply(seq_len(total_i), function(i){
    print(round(i/total_i*100, 2))
    sample_t <- data_ori$Sample[i]
    pos_t <- data_ori$Position[i]
    seq_t <- hk_seqs_1623[names(hk_seqs_1623) == sample_t]
    con_base_t <- subseq(seq_t, pos_t, pos_t)
    check <- con_base_t == data_ori$`Alt base`[i]
    if(check){return("Major")}else{return("Minor")} # The Major variants refer to the SNVs that presented on the consensus sequences. 
})

write_csv(data_ori, "../../Bottleneck_and_SNVs/data_ori_all_variants_001.csv")
data_ori <- data_ori %>% filter(as.numeric(Alle_freq) >= 0.03) # supported by reapeated samples 
write_csv(data_ori, "../../Bottleneck_and_SNVs/data_ori_all_variants.csv")

# data_ori_minor <- data_ori %>% filter(as.numeric(Alle_freq) < 0.95)
# write_csv(data_ori_minor, "../../Bottleneck_and_SNVs/data_ori_minor_variants.csv")


