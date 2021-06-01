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
library(ggpubr)

# mc.cores=12
mc.cores=48

identify_sig <- function(data, group, value){
    tmp <- combn(unique(as.character(data[[group]])), 2)
    tmp <- lapply(seq_len(ncol(tmp)), function(j){
        # print(j)
        pair_t <- tmp[,j]
        value1 <- data[[value]][data[[group]] == pair_t[1]]
        value2 <- data[[value]][data[[group]] == pair_t[2]]
        check <- wilcox.test(value1, value2)
        if(check$p.value<=0.05){
            rst_out <- pair_t
            names(rst_out) <- rep(round(check$p.value, 3), length(pair_t))
            return(rst_out)
        } else {
            return(NA)
        }
    })
    return(tmp[!is.na(tmp)])
}

ref_seq <- readDNAStringSet("../../2020-09-01_COVID_NGS_pipeline/NGS_data_input/reference.fasta")
hk_seqs <- readDNAStringSet("../results/hk_case.fasta")
hk_seqs_1623 <- readDNAStringSet("../results/hk_seq_1623.fasta")
seqed_meta <- read_csv("../data/cov_and_lin.csv") 
id_case_df <- read_csv("../data/labID_caseID_merged.csv", col_types = cols(.default = "c"))
data_metadata <- readxl::read_excel("../data/Copy of 2019ncov-global-surveillance-linelist (2021.01.26).xlsx", col_types = "text")
data_metadata <- data_metadata %>% select(`HK case no.`:`onset to report (days)`)
data_metadata <- data_metadata[!is.na(data_metadata$`HK case no.`),]
data_metadata$case_id <- seq_len(nrow(data_metadata))
data_metadata$club_sim <- gsub("[^[:alnum:]]", "_", data_metadata$`Bar/pub/club`)
data_metadata <- data_metadata[!duplicated(data_metadata$`HK case no.`),]
metadata <- data_metadata
id_case_df_hk <- id_case_df %>% filter(case_ID %in% names(hk_seqs))
data_hk_study_samples <- read_csv("../results/samples_hk_study.csv")
df_sample_type <- read_csv("../results/wave_type_rp.csv")
data_duplicated_samples <- read_csv("../results/duplicated_cases.csv")
data_duplicated_samples_keep <- data_duplicated_samples %>% filter(!(same_sample & diff_in_consensus>0)) 
data_duplicated_samples_keep <- data_duplicated_samples_keep %>% filter(`Genome coverage`*29903>=27000)
data_duplicated_samples_keep <- data_duplicated_samples_keep[data_duplicated_samples_keep$case_id %in% data_duplicated_samples_keep$case_id[duplicated(data_duplicated_samples_keep$case_id)],]
data_duplicated_samples_keep <- data_duplicated_samples_keep %>% filter(!case_id == "7409")
df_pair <- read_csv("../results/identified_transmission_pairs.csv")
data_ori <- read_csv("../results/data_ori_all_variants.csv")

sample_should_have <- c(data_hk_study_samples$sample, data_duplicated_samples_keep$`Lab No.`)

# data duplicated
data <- data_ori[data_ori$Sample %in% data_duplicated_samples_keep$Sample, ]
quantile(data$depth, seq(0,1, 0.01))

## overalapped mutation for transmission pairs
data <- data_ori
df_pair_day1 <- df_pair %>% filter(min_date_diff == 1)
df_pair_plot <- lapply(seq_len(nrow(df_pair_day1)), function(i){
	case_id_t <- df_pair_day1$case_id[i]
    case_id_t <- strsplit(case_id_t, " ")[[1]]
    sample_t <- id_case_df_hk %>% filter(case_ID %in% case_id_t) %>% .$lab_ID
    		
	# data_filter <- data %>% filter(Alle_freq < 0.95)
    data_filter <- data_ori
	data_snvs_tmp <- data_filter %>% filter(case_id %in% case_id_t)
    data_snvs_tmp_major <- data_snvs_tmp %>% filter(iSNVs == "Major")
    data_snvs_tmp_isnvs <- data_snvs_tmp %>% filter(iSNVs == "Minor")

	mut_common_major <- (table(data_snvs_tmp_major$mutation) == 2)
	mut_common_major <- unique(names(mut_common_major)[mut_common_major])
	mut_self_major <- (table(data_snvs_tmp_major$mutation) == 1)
	mut_self_major <- unique(names(mut_self_major)[mut_self_major])
    
    mut_common_isnvs <- (table(data_snvs_tmp_isnvs$mutation) == 2)
	mut_common_isnvs <- unique(names(mut_common_isnvs)[mut_common_isnvs])
	mut_self_isnvs <- (table(data_snvs_tmp_isnvs$mutation) == 1)
	mut_self_isnvs <- unique(names(mut_self_isnvs)[mut_self_isnvs])

    return(tibble(case_id  = df_pair$case_id[i], `Shared major SNVs` = length(mut_common_major), `Unique major SNVs` = length(mut_self_major), `Shared minor SNVs` = length(mut_common_isnvs), `Unique minor SNVs` = length(mut_self_isnvs)))
})
df_pair_plot <- bind_rows(df_pair_plot)
df_pair_plot <- left_join(df_pair_day1, df_pair_plot)
df_pair_plot$pair_sim <- gsub("Cluster_", "", df_pair_plot$pair)

colors <- c("Shared major" = "#ca0020", "Shared minor" = "#0571b0", "Unique major" = "#f4a582", "Unique minor" = "#92c5de")
(p1 <- df_pair_plot %>% pivot_longer(contains("Shared") | contains("Unique")) %>% mutate(name = gsub(" SNVs", "", name)) %>% ggplot(aes(x = pair_sim, y = value))+
    geom_col(aes(fill = name), position = "fill", width = 0.618) +
    # scale_fill_uchicago(name = "Type")+
    scale_fill_manual(name = "SNV type", values = colors, labels = c("Major (shared)", "Minor (shared)", "Major (unique)", "Minor (unique)"))+
    xlab("Transmission pair")+
    ylab("Proportion")+
    scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    theme_classic())

## Alle frequency between donor and recipient
df_pair_day1_full <- read_csv("../results/bottleneck_1day.csv")
df_pair_day1_full <- df_pair_day1_full %>% filter(var_calling_threshold==0.03)
df_pair_day1_full_plot <- df_pair_day1_full %>% arrange(pair) %>% left_join(df_pair_day1)

df_pair_day1_full_plot$pair_sim <- gsub("Cluster_", "", df_pair_day1_full_plot$pair)

colors <- c("#ffffcc", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#0c2c84")
(p2 <- df_pair_day1_full_plot %>% ggplot(aes(x = pair_sim, y = bottleneck_size, fill = factor(date_diff)))+
    geom_errorbar(aes(ymin=CI_index_lower, ymax=CI_index_upper), color = "black", width = 0.2)+
    geom_point(size = 3, shape = 21, alpha = 0.9)+
    ylim(0,220)+
    scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    scale_fill_manual(name = "Lag (days)", values = colors)+
    xlab("Transmission pair")+
    ylab("Bottlenect size")+
    theme_classic())

(p3 <- ggarrange(p1, p2, ncol = 1, align = "v"))
ggsave(paste0("../results/transmission_pair_mut_bottleneck.pdf"), plot = p3, width = 8, height = 6)


# Jaccard distance of all sample pairs
data <- data_ori[data_ori$Sample %in% c(sample_should_have, data_duplicated_samples_keep$Sample), ]
df_caseid_labid <- data %>% select(Sample, case_id) %>% unique()
combn_samples <- combn(df_caseid_labid$Sample, 2)
dim(combn_samples)

df_jd_all <- mclapply(seq_len(ncol(combn_samples)), function(j){
    print(j)
    combn_tmp <- combn_samples[,j]
    data_tmp <- data %>% filter(Sample %in% combn_tmp)
    if(nrow(data_tmp)<1){
        return(NA)
    } else{
        data_tmp_1 <- data_tmp %>% filter(Sample == combn_tmp[1])
        data_tmp_2 <- data_tmp %>% filter(Sample == combn_tmp[2])

        data_tmp_1_major <- data_tmp_1 %>% filter(iSNVs == "Major")
        data_tmp_1_minor <- data_tmp_1 %>% filter(iSNVs == "Minor")
        data_tmp_2_major <- data_tmp_2 %>% filter(iSNVs == "Major")
        data_tmp_2_minor <- data_tmp_2 %>% filter(iSNVs == "Minor")

        mutation_union_n <- length(unique(union(data_tmp_1$mutation, data_tmp_2$mutation)))
        mutation_instersect_n <- length(intersect(data_tmp_1$mutation, data_tmp_2$mutation))
        jd_all <- 1-mutation_instersect_n/mutation_union_n

        mutation_union_major <- length(unique(union(data_tmp_1_major$mutation, data_tmp_2_major$mutation)))
        mutation_instersect_major <- length(intersect(data_tmp_1_major$mutation, data_tmp_2_major$mutation))
        jd_major <- 1-mutation_instersect_major/mutation_union_major

        mutation_union_minor <- length(unique(union(data_tmp_1_minor$mutation, data_tmp_2_minor$mutation)))
        mutation_instersect_minor <- length(intersect(data_tmp_1_minor$mutation, data_tmp_2_minor$mutation))
        jd_minor <- 1-mutation_instersect_minor/mutation_union_minor

        return(tibble(Sample_1 = combn_tmp[1], Sample_2 = combn_tmp[2], case_id_1 = df_caseid_labid$case_id[df_caseid_labid$Sample == combn_tmp[1]], case_id_2 = df_caseid_labid$case_id[df_caseid_labid$Sample == combn_tmp[2]], mutation_union_all = mutation_union_n, mutation_instersect_all = mutation_instersect_n, jd_all = jd_all, mutation_union_major = mutation_union_major, mutation_instersect_major = mutation_instersect_major, jd_major = jd_major, mutation_union_minor = mutation_union_minor, mutation_instersect_minor = mutation_instersect_minor, jd_minor = jd_minor))
    }
}, mc.cores = 24)
df_jd_all <- df_jd_all[!is.na(df_jd_all)]
df_jd_all <- bind_rows(df_jd_all)
write_csv(df_jd_all, "../results/df_jd_all.csv")
rm(combn_samples)
df_jd_all <- read_csv("../results/df_jd_all.csv")

position_snps <- unique(data$Position)

#### check the intra-, inter- clsuter jd
data_cluster <- read_csv("../results/data_cluster.csv")

df_jd_all$case_id_1 <- factor(df_jd_all$case_id_1)
df_jd_all$case_id_2 <- factor(df_jd_all$case_id_2)

wave_t <- sapply(levels(df_jd_all$case_id_1), function(x){
    tmp <- df_sample_type$Type[df_sample_type$case_id==x][1]
    if(length(tmp)==0){
        return(NA)
    } else {
        return(tmp)
    }
})
df_jd_all$wave_1 <- wave_t[df_jd_all$case_id_1]

wave_t <- sapply(levels(df_jd_all$case_id_2), function(x){
    tmp <- df_sample_type$Type[df_sample_type$case_id==x][1]
    if(length(tmp)==0){
        return(NA)
    } else {
        return(tmp)
    }
})
df_jd_all$wave_2 <- wave_t[df_jd_all$case_id_2]

df_jd_all$samewave <- df_jd_all$wave_1 == df_jd_all$wave_2 
df_jd_all$samewave <- ifelse(df_jd_all$samewave, "Within same wave", "Not in same wave")
df_jd_all$samewave[is.na(df_jd_all$samewave)] <- "Not in same wave"

cluster_1 <- sapply(levels(df_jd_all$case_id_1), function(x){
    tmp <- data_cluster$metaCluster[data_cluster$case_id==x][1]
    if(length(tmp)==0){
        return(NA)
    } else {
        return(tmp)
    }
})
df_jd_all$cluster_1 <- cluster_1[df_jd_all$case_id_1]

cluster_2 <- sapply(levels(df_jd_all$case_id_2), function(x){
    tmp <- data_cluster$metaCluster[data_cluster$case_id==x][1]
    if(length(tmp)==0){
        return(NA)
    } else {
        return(tmp)
    }
})
df_jd_all$cluster_2 <- cluster_2[df_jd_all$case_id_2]

df_jd_all$metaCluster_check <- df_jd_all$cluster_1 == df_jd_all$cluster_2 
df_jd_all$metaCluster_check <- ifelse(df_jd_all$metaCluster_check, "Within same cluster", "Not in same cluster")
df_jd_all$metaCluster_check[is.na(df_jd_all$metaCluster_check)] <- "Not in same cluster"

minorCluster <- mclapply(seq_len(nrow(df_jd_all)), function(i){
	print(i)
    tmp1 <- data_cluster$value[data_cluster$case_id == df_jd_all$case_id_1[i]]
    tmp2 <- data_cluster$value[data_cluster$case_id == df_jd_all$case_id_2[i]]
    
    if(length(tmp1)==0 | length(tmp2)==0){
        return(NA)
    }

	tmp_t <- intersect(tmp1, tmp2)

    if(length(tmp_t)==0){
        return("Different")
    } else{
        return(tmp_t[1])
    }
}, mc.cores = 24)

df_jd_all$minorCluster <- unlist(minorCluster)
df_jd_all$minorCluster_check <- df_jd_all$minorCluster
df_jd_all$minorCluster_check[df_jd_all$minorCluster == "Different"] <- "Without direct link"
df_jd_all$minorCluster_check[is.na(df_jd_all$minorCluster_check)] <- "Without direct link"
df_jd_all$minorCluster_check[df_jd_all$minorCluster_check != "Without direct link"] <- "With direct link"

df_jd_all$meta_minor_cluster <- df_jd_all$minorCluster_check
df_jd_all$meta_minor_cluster[df_jd_all$minorCluster_check == "Without direct link" & df_jd_all$metaCluster_check == "Not in same cluster"] <- "Not in same cluster"
df_jd_all$meta_minor_cluster[df_jd_all$meta_minor_cluster == "Without direct link"] <- "Same cluster without direct link"

df_jd_all$same_sample_check <- as.character(df_jd_all$case_id_1) == as.character(df_jd_all$case_id_2)
df_jd_all$id <- paste(df_jd_all$Sample_1, df_jd_all$Sample_2)

write_csv(df_jd_all, "../results/df_jd_all.csv")
df_jd_all <- read_csv("../results/df_jd_all.csv")

df_jd_all$meta_minor_cluster_same_patient <- df_jd_all$meta_minor_cluster
df_jd_all$meta_minor_cluster_same_patient[as.character(df_jd_all$case_id_1) == as.character(df_jd_all$case_id_2)] <- "Same patient"

check1 <- paste(df_jd_all$case_id_1, df_jd_all$case_id_2) %in% df_pair_day1$case_id
check2 <- paste(df_jd_all$case_id_2, df_jd_all$case_id_1) %in% df_pair_day1$case_id
df_jd_all[check1 | check2,]
df_jd_all$meta_minor_cluster_same_patient[check1 | check2] <- "Identified transmission pairs"
df_jd_all$meta_minor_cluster_same_patient <- factor(df_jd_all$meta_minor_cluster_same_patient, levels = c("Not in same cluster", "Same cluster without direct link", "With direct link", "Identified transmission pairs", "Same patient"), labels = c("Not in same cluster", "Same cluster (without direct link)", "Same cluster (with direct link)", "Identified transmission pairs", "Same patient"))

df_jd_all$meta_minor_cluster_same_patient_sim <- as.character(df_jd_all$meta_minor_cluster_same_patient)
df_jd_all$meta_minor_cluster_same_patient_sim[grepl("^Same cluster", df_jd_all$meta_minor_cluster_same_patient_sim)] <- "Same cluster"
df_jd_all$meta_minor_cluster_same_patient_sim <- factor(df_jd_all$meta_minor_cluster_same_patient_sim, levels = c("Not in same cluster", "Same cluster", "Identified transmission pairs", "Same patient"))

df_jd_all %>% group_by(meta_minor_cluster_same_patient_sim) %>% summarise(mean(jd_major, na.rm = T))
tmp1 <- identify_sig(df_jd_all, "meta_minor_cluster_same_patient_sim", "jd_major")
p4_1 <- ggplot(df_jd_all, aes(x = meta_minor_cluster_same_patient_sim, y = jd_major))+
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.8) + 
    geom_signif(comparisons = tmp1, map_signif_level = TRUE, textsize = 6, step_increase = 0.1)+
    scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    ylab("Jaccard distance of major SNVs") + xlab("Sample pairs")+
    theme_classic()

df_jd_all %>% group_by(meta_minor_cluster_same_patient_sim) %>% summarise(mean(jd_minor, na.rm = T))
tmp2 <- identify_sig(df_jd_all, "meta_minor_cluster_same_patient_sim", "jd_minor")
p4_2 <- ggplot(df_jd_all, aes(x = meta_minor_cluster_same_patient_sim, y = jd_minor))+
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.8) + 
    geom_signif(comparisons = tmp2, map_signif_level = TRUE, textsize = 6, step_increase = 0.1)+
    scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    ylab("Jaccard distance of minor SNVs") + xlab("Sample pairs")+
    theme_classic()

p4 <- ggarrange(p4_1, p4_2, nrow = 1, align = "hv")
ggsave("../results/jaccard_dist_supple.pdf", width = 12/1.414, height = 12/2, plot = p4)

# rain cloud plot
library(raincloudplots)

data_2x2_mod <- function (array_1, array_2, array_3, array_4, array_5, array_6, array_7, array_8, labels, jit_distance = 0.05, jit_seed = 2021, spread_x_ticks = TRUE) {
    n <- max(length(array_1), length(array_2), length(array_3), 
        length(array_4), length(array_5), length(array_6), length(array_7), length(array_8))
    length(array_1) <- n
    length(array_2) <- n
    length(array_3) <- n
    length(array_4) <- n
    length(array_5) <- n
    length(array_6) <- n
    length(array_7) <- n
    length(array_8) <- n
    
    if (spread_x_ticks == TRUE){
        data_2x2 <- data.frame(
            y_axis = c(array_1, array_2, array_3, array_4, array_5, array_6, array_7, array_8), 
            x_axis = rep(c(1, 1.01, 2, 2.01, 3, 3.01, 4, 4.01), each = n),
            id = factor(rep(1:n, 2)),
            group = rep(c(labels[1], labels[2]), 
            each = n))
    } 
        
    data_2x2$jit <- jitter(data_2x2$x_axis, amount = jit_distance)
    if (any(is.na(data_2x2))) 
        data_2x2 <- stats::na.omit(data_2x2)
    return(data_2x2)
}

raincloud_2x3_repmes <- function(data_2x2,
                                 colors = (c('#ca0020', '#0571b0', '#ca0020', '#0571b0','#ca0020', '#0571b0','#ca0020', '#0571b0')),
                                 fills = (c('#ca0020', '#0571b0','#ca0020', '#0571b0','#ca0020', '#0571b0','#ca0020', '#0571b0')),
                                 size = 1.5,
                                 alpha = .8,
                                 alpha_point = .1) {
    figure_2x3 <- ggplot(data_2x2)
    #Add geom_() objects
    # for(i in 1:4){
    #     data_tmp <- data_2x2 %>% dplyr::filter(x_axis == i) %>% sample_n(min(2000, n()))
    #     alpha_tmp <- ifelse(nrow(data_tmp)>50, alpha_point, 0.6)
    #     figure_2x3 <- figure_2x3 + geom_point(data = data_tmp, aes(x = jit, y = y_axis), color = colors[i*2-1], fill = fills[i*2-1], size = size, alpha = alpha_tmp)

    #     data_tmp <- data_2x2 %>% dplyr::filter(x_axis == i+0.01) %>% sample_n(min(2000, n()))
    #     alpha_tmp <- ifelse(nrow(data_tmp)>50, alpha_point, 0.6)
    #     figure_2x3 <- figure_2x3 + geom_point(data = data_tmp, aes(x = jit, y = y_axis), color = colors[i*2], fill = fills[i*2], size = size, alpha = alpha_tmp)
    # }

    for(i in 1:4){
        figure_2x3 <- figure_2x3 + geom_half_boxplot(
        data = data_2x2 %>% dplyr::filter(x_axis==i), aes(x=x_axis, y = y_axis, fill = fills[i*2-1], color = colors[i*2-1]), position = position_nudge(x = .2),
        side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, alpha = alpha)
       
        figure_2x3 <- figure_2x3 + geom_half_boxplot(
        data = data_2x2 %>% dplyr::filter(x_axis==i+0.01), aes(x=x_axis, y = y_axis, fill = fills[i*2], color = colors[i*2]), position = position_nudge(x = .2),
        side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, alpha = alpha)

    }

    for(i in 1:4){
        figure_2x3 <- figure_2x3 + 
        geom_half_violin(
        data = data_2x2 %>% dplyr::filter(x_axis==i),aes(x = x_axis, y = y_axis), color = colors[i*2-1], fill = fills[i*2-1], position = position_nudge(x = .45), side = "r", alpha = alpha) 
           
        figure_2x3 <- figure_2x3 + geom_half_violin(
        data = data_2x2 %>% dplyr::filter(x_axis==i+0.01),aes(x = x_axis, y = y_axis), color = colors[i*2], fill = fills[i*2], position = position_nudge(x = .45), side = "r", alpha = alpha) 

    }

    return(figure_2x3)

}

levels(df_jd_all$meta_minor_cluster_same_patient_sim)

df_2_4 <- data_2x2_mod(
    array_1 = df_jd_all$jd_major[df_jd_all$meta_minor_cluster_same_patient_sim=="Not in same cluster"],
    array_2 = df_jd_all$jd_minor[df_jd_all$meta_minor_cluster_same_patient_sim=="Not in same cluster"],
    array_3 = df_jd_all$jd_major[df_jd_all$meta_minor_cluster_same_patient_sim=="Same cluster"],
    array_4 = df_jd_all$jd_minor[df_jd_all$meta_minor_cluster_same_patient_sim=="Same cluster"],
    array_5 = df_jd_all$jd_major[df_jd_all$meta_minor_cluster_same_patient_sim=="Identified transmission pairs"],
    array_6 = df_jd_all$jd_minor[df_jd_all$meta_minor_cluster_same_patient_sim=="Identified transmission pairs"],
    array_7 = df_jd_all$jd_major[df_jd_all$meta_minor_cluster_same_patient_sim=="Same patient"],
    array_8 = df_jd_all$jd_minor[df_jd_all$meta_minor_cluster_same_patient_sim=="Same patient"],
    labels = c('Major SNVs','Minor SNVs')
    )

p5_0 <- raincloud_2x3_repmes(df_2_4, size = 1)

(p5_1 <- p5_0 + scale_x_continuous(
    breaks=c(1.4,2.4,3.4,4.4),
    labels=c("Not in same cluster", "Same cluster", "Identified transmission pairs", "Same patient"),
    # guide = guide_axis(n.dodge = 2),
    limits=c(1,5)
    ) +
  xlab("Sample pairs") + 
  ylab("Jaccard distance") +
#   coord_flip()+
  theme_classic())


tmp1 <- identify_sig(df_jd_all, "meta_minor_cluster_same_patient_sim", "jd_major")
tmp2 <- identify_sig(df_jd_all, "meta_minor_cluster_same_patient_sim", "jd_minor")
tmp1 # significant pair, major SNVs
tmp2 # significant pair, minor SNVs

p5_2 <- p5_1
cur_y = 1.05
sapply(tmp1, function(x){
    idx <- sapply(x, function(y){
        which(y == levels(df_jd_all$meta_minor_cluster_same_patient_sim))
    })
    idx <- idx + 0.4
    label_t <- names(idx)[1]
    if(label_t == "0"){label_t <- "***"}
    p5_2 <<- p5_2 +
        annotate("segment", 
            x=c(idx[1],idx[1],idx[2]),
            xend=c(idx[1],idx[2],idx[2]),
            y= c(cur_y-0.01,cur_y,cur_y),
            yend=c(cur_y,cur_y,cur_y-0.01),
            color = "#ca0020")+
        annotate("text", 
            x=mean(idx),
            y= cur_y + 0.015,
            label = label_t,
            size = 2.5,
            color = "#ca0020")
    cur_y <<- cur_y + 0.03
})

sapply(tmp2, function(x){
    idx <- sapply(x, function(y){
        which(y == levels(df_jd_all$meta_minor_cluster_same_patient_sim))
    })
    idx <- idx + 0.4
    label_t <- names(idx)[1]
    if(label_t == "0"){label_t <- "***"}
    p5_2 <<- p5_2 +
        annotate("segment", 
            x=c(idx[1],idx[1],idx[2]),
            xend=c(idx[1],idx[2],idx[2]),
            y= c(cur_y-0.01,cur_y,cur_y),
            yend=c(cur_y,cur_y,cur_y-0.01),
            color = "#0571b0")+
        annotate("text", 
            x=mean(idx),
            y= cur_y + 0.015,
            label = label_t,
            size = 2.5,
            color = "#0571b0")
    cur_y <<- cur_y + 0.03
})

(p5_3 <- p5_2 + 
    scale_fill_identity(name = "SNV type", breaks = c("#ca0020", "#0571b0"),labels = c("Major SNVs", "Minor SNVs"), guide = "legend")+
    scale_color_identity(name = "SNV type", breaks = c("#ca0020", "#0571b0"),labels = c("Major SNVs", "Minor SNVs"), guide = "legend"))

p_out <- ggarrange(p3, p5_3, 
    # heights = c(0.55, 0.45),
    labels = c("A", "B"),
    ncol = 1)
ggsave("../results/main_figure_bottleneck_snvs.pdf", width = 12/1.414, height = 12, plot = p_out)
