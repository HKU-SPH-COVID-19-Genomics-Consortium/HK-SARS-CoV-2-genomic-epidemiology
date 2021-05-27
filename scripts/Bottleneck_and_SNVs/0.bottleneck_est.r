# This code requires use of the R package rmutil
library(rmutil)
# library(argparse)
# library(mdatools)
library(tidyverse)
library(readxl)
library(MASS)
library(parallel)
library(Biostrings)
library(ggsci)
library(scales)

# mc.cores=12
mc.cores=48

##########################################################################

generate_log_likelihood_function_approx <- function(donor_freqs_observed, recipient_total_reads, recipient_var_reads_observed, Nb_min, Nb_max, var_calling_threshold, confidence_level, n_variants) {
  num_NB_values <- Nb_max - Nb_min + 1
  log_likelihood_function <- matrix(0, Nb_max)
  log_likelihood_matrix_t <- matrix(0, n_variants, num_NB_values)
  sapply(seq_len(n_variants), function(i){
    rst <- mclapply(seq_len(num_NB_values), function(j){
      Nb_val <- (j - 1 + Nb_min)
      nu_donor <- donor_freqs_observed[i, 1]
      variant_reads <- recipient_var_reads_observed[i, 1]
      total_reads <- recipient_total_reads[i, 1]
      nu_recipient <- variant_reads/total_reads
      tmp <- 0
      if(nu_recipient == 1){nu_recipient <- 0.999}
      if (nu_recipient >= var_calling_threshold) { # implement variant calling threshold
        for (k in 0:Nb_val) {
          # print(tmp)
          tmp <- tmp + (dbeta(nu_recipient, k, (Nb_val - k)) * dbinom(k, size = Nb_val, prob = nu_donor))
        }
      }
      if (nu_recipient < var_calling_threshold) { # implement variant calling threshold
        for (k in 0:Nb_val) {
          tmp <- tmp + (pbeta(var_calling_threshold, k, (Nb_val - k)) * dbinom(k, size = Nb_val, prob = nu_donor))
        }
      }
      return(log(tmp))
    }, mc.cores = mc.cores)
    # Now we sum over log likelihoods of the variants at different loci to get the total log likelihood for each value of Nb
    log_likelihood_function[Nb_min:length(log_likelihood_function)] <<- log_likelihood_function[Nb_min:length(log_likelihood_function)] + unlist(rst)
    log_likelihood_matrix_t[i,] <<- unlist(rst)
    print(paste0(i, " of ", n_variants))
  })
  return(log_likelihood_function)
}

restrict_log_likelihood <- function(log_likelihood_function, Nb_min, Nb_max) # restricts log likelihood to the interval of interst
{
  for (h in 1:(Nb_min)) {
    if (h < Nb_min) {
      log_likelihood_function[h] <- -999999999
    } # kludge for ensuring that these values less than Nb_min don't interfere with our search for the max of log likelihood in the interval of Nb_min to Nb_max
  }

  return(log_likelihood_function)
  # print(erfinv(percent_confidence_interval)*sqrt(2))
}

return_bottleneck_size <- function(log_likelihood_function, Nb_min, Nb_max) {
  max_log_likelihood <- which(log_likelihood_function == max(log_likelihood_function))

  return(max_log_likelihood)
}

erfinv <- function(x) qnorm((1 + x) / 2) / sqrt(2)

return_CI_lower <- function(log_likelihood_function, Nb_min, Nb_max) ## returns lower bound of confidence interval
{
  max_log_likelihood <- which(log_likelihood_function == max(log_likelihood_function)) ## This is the point on the x-axis (bottleneck size) at which log likelihood is maximized
  max_val <- max(log_likelihood_function) ## This is the maximum value of the log likelihood function, found when the index is our bottleneck estimate
  CI_height <- max_val - erfinv(confidence_level) * sqrt(2) # This value (  height on y axis) determines the confidence intervals using the likelihood ratio test


  CI_index_lower <- Nb_min
  CI_index_upper <- max_log_likelihood
  for (h in 1:Nb_min) {
    if (h < Nb_min) {
      log_likelihood_function[h] <- NA
    } #  Removing parameter values less than Nb_min from plot
  }
  ## above loop just enforces our minimum bottleneck cutoff
  for (h in Nb_min:max_log_likelihood) {
    test1 <- (log_likelihood_function[CI_index_lower] - CI_height) * (log_likelihood_function[CI_index_lower] - CI_height)
    test2 <- (log_likelihood_function[h] - CI_height) * (log_likelihood_function[h] - CI_height)
    if (test2 < test1) {
        # print(h)
      CI_index_lower <- h
    }
  }
  if ((log_likelihood_function[CI_index_lower] - CI_height) > 0) {
    CI_index_lower <- CI_index_lower - 1
  }
  # above loops use likelihood ratio test to find lower confidence interval
  for (h in max_log_likelihood:Nb_max)
  {
    test1 <- (log_likelihood_function[CI_index_upper] - CI_height) * (log_likelihood_function[CI_index_upper] - CI_height)
    test2 <- (log_likelihood_function[h] - CI_height) * (log_likelihood_function[h] - CI_height)
    if (test2 < test1) {
      CI_index_upper <- h
    }
  }
  if ((log_likelihood_function[CI_index_upper] - CI_height) > 0) {
    CI_index_upper <- CI_index_upper + 1
  }

  return(CI_index_lower)
}


return_CI_upper <- function(log_likelihood_function, Nb_min, Nb_max) ## returns upper bound of confidence interval
{
  max_log_likelihood <- which(log_likelihood_function == max(log_likelihood_function)) ## This is the point on the x-axis (bottleneck size) at which log likelihood is maximized
  max_val <- max(log_likelihood_function) ## This is the maximum value of the log likelihood function, found when the index is our bottleneck estimate
  CI_height <- max_val - erfinv(confidence_level) * sqrt(2) # This value (  height on y axis) determines the confidence intervals using the likelihood ratio test

  CI_index_lower <- Nb_min
  CI_index_upper <- max_log_likelihood
  for (h in 1:Nb_min) {
    if (h < Nb_min) {
      log_likelihood_function[h] <- NA
    } #  Removing parameter values less than Nb_min from plot
  }
  ## above loop just enforces our minimum bottleneck cutoff
  for (h in Nb_min:max_log_likelihood) {
    test1 <- (log_likelihood_function[CI_index_lower] - CI_height) * (log_likelihood_function[CI_index_lower] - CI_height)
    test2 <- (log_likelihood_function[h] - CI_height) * (log_likelihood_function[h] - CI_height)
    if (test2 < test1) {
      CI_index_lower <- h
    }
  }
  if ((log_likelihood_function[CI_index_lower] - CI_height) > 0) {
    CI_index_lower <- CI_index_lower - 1
  }
  # above loops use likelihood ratio test to find lower confidence interval
  for (h in max_log_likelihood:Nb_max)
  {
    test1 <- (log_likelihood_function[CI_index_upper] - CI_height) * (log_likelihood_function[CI_index_upper] - CI_height)
    test2 <- (log_likelihood_function[h] - CI_height) * (log_likelihood_function[h] - CI_height)
    if (test2 < test1) {
      CI_index_upper <- h
    }
  }
  if ((log_likelihood_function[CI_index_upper] - CI_height) > 0) {
    CI_index_upper <- CI_index_upper + 1
  }

  return(CI_index_upper)
}

### identify transimission pairs
data_alignment <- readDNAStringSet("../../data/seq_local_analysis/hk_case_primer_masked_2021-05-08.afa")
data_hk_study_samples <- read_csv("../../data/Bottleneck_and_SNVs/samples_hk_study.csv")
df_raw_data <- read_csv("../../data/Bottleneck_and_SNVs/cov_and_lin.csv")
data_metadata <- readxl::read_excel("../../data/Copy of 2019ncov-global-surveillance-linelist (2021.01.26).xlsx", col_types = "text")
data_metadata <- data_metadata[!is.na(data_metadata$`HK case no.`),]
data_metadata$case_id <- seq_len(nrow(data_metadata))
data_metadata$club_sim <- gsub("[^[:alnum:]]", "_", data_metadata$`Bar/pub/club`)
data_metadata <- data_metadata[!duplicated(data_metadata$`HK case no.`),]
data <- data_metadata
data_cluster <- data %>% dplyr::select(case_id, Cluster:`Quaternary generation 1`) %>% filter(Cluster == "Y") %>% pivot_longer(`Cluster setting (primary 1)`:`Quaternary generation 1`) %>% filter(!is.na(value))

#### filtering criteria: 
#### 1. the receipient case was included in only one minor cluster, 
#### 2. the minor cluster only included two cases,
#### 3. the onset date of the recipient case is at aleast five days after the onset date of donor case.
#### 4. we had sequenced both the donor and recipient cases
metaCluster_two_cases <- names(table(data$metaCluster)[table(data$metaCluster) == 2])
identify_transmission_pairs <- function(data_cluster=data_cluster, data=data, min_date_diff = 5){
  tmp <- table(data_cluster$value)
  trans_pairs <- names(tmp[tmp == 2]) # two cases
  case_id_t <- c()
  trans_pairs_filter <- sapply(trans_pairs, function(x){
      print(x)
      df_tmp <- data_cluster %>% filter(value == x)
      # df_tmp <- data_cluster %>% filter(case_id %in% df_tmp$case_id)
      # if(nrow(df_tmp)!=2){  # two cases without other links
      #   return(NA)
      # }
      data_tmp <- data %>% filter(case_id %in% df_tmp$case_id)
      if(!all(data_tmp$metaCluster %in% metaCluster_two_cases)){ # two cases without other links
        return(NA)
      }
      date_tmp <- lubridate::ymd(data_tmp$`Onset date`)  
      date_diff <- as.numeric(date_tmp[which.max(date_tmp)] - date_tmp[-which.max(date_tmp)])
      if(length(date_diff)==0){
          return(NA)
      }
      if(is.na(date_diff)){
          return(NA)
      }
      check1 <- date_diff >= min_date_diff # >= 5 days
      if(check1){
          check2 <- sum(data_cluster$case_id == data_tmp$case_id[which.max(date_tmp)]) < 2 # recipient only in one minor cluster
          if(check2){
              check3 <- all(data_tmp$case_id %in% df_raw_data$`HK case no.`) & all(data_tmp$case_id %in% names(data_alignment))
              if(check3){
                  case_id_t <<- c(case_id_t, list(data_tmp$case_id))
                  return(date_diff)
              } else {
                  return(NA)
              }
          } else {
              return(NA)
          }
      } else {
          return(NA)
      }
  })
  trans_pairs <- trans_pairs[!is.na(trans_pairs_filter)] 
  date_diff <- trans_pairs_filter[!is.na(trans_pairs_filter)]
  case_id_t  <- sapply(case_id_t, function(x){paste(x, collapse = " ")})
  return(tibble(pair = trans_pairs, min_date_diff = min_date_diff, case_id = case_id_t, date_diff = date_diff)) 
}
pair_1 <- identify_transmission_pairs(data_cluster=data_cluster, data=data, min_date_diff = 1)
pair_2 <- identify_transmission_pairs(data_cluster=data_cluster, data=data, min_date_diff = 2)
pair_3 <- identify_transmission_pairs(data_cluster=data_cluster, data=data, min_date_diff = 3)
pair_4 <- identify_transmission_pairs(data_cluster=data_cluster, data=data, min_date_diff = 4)
pair_5 <- identify_transmission_pairs(data_cluster=data_cluster, data=data, min_date_diff = 5)
df_pair <- bind_rows(pair_1, pair_2, pair_3, pair_4, pair_5)
# tmp <- read_csv("../results/identified_transmission_pairs.csv")
# tmp$pair %in% df_pair$pair
# df_pair[!df_pair$pair %in% tmp$pair,]
write_csv(df_pair, "../../results/Bottleneck_and_SNVs/identified_transmission_pairs.csv")
df_pair <- read_csv("../../results/Bottleneck_and_SNVs/identified_transmission_pairs.csv")

df_masked_sites <- read_csv("../../data/Bottleneck_and_SNVs/primer_masked_sites.csv")
masked_sites <- c(df_masked_sites$Position, 1:100, (29903-99):29903)
df_bottleneck_all <- tibble()

for(day in 1){
  var_calling_threshold_all <- c(0.03, 0.06)
  # var_calling_threshold_all <- c(0.03)
  Nb_min <- 1
  Nb_max <- 1000
  confidence_level <- 0.95
  trans_pairs <- df_pair$pair[df_pair$min_date_diff == day]

  df_trans_pair_reads <- tibble()
  tmp <- rep(trans_pairs, length(var_calling_threshold_all))
  df <- tibble(pair = tmp,
      var_calling_threshold = rep(var_calling_threshold_all, each = length(trans_pairs)),
      Donor = rep(NA, length(tmp)), 
      Recipient = rep(NA, length(tmp)), 
      bottleneck_size = rep(NA, length(tmp)), 
      CI_index_lower = rep(NA, length(tmp)), 
      CI_index_upper = rep(NA, length(tmp)))
  log_likelihood_matrix <- matrix(0, length(tmp),  Nb_max - Nb_min + 1)

  sapply(seq_len(nrow(df)), function(i){
      id <- df$pair[i]
      var_calling_threshold <- df$var_calling_threshold[i]
      print("-------------------------")
      print(id)
      print(var_calling_threshold)
      print("-------------------------")
      
      df_tmp <- data_cluster %>% filter(value == id)
      data_tmp <- data %>% filter(case_id %in% df_tmp$case_id)
      date_tmp <- lubridate::ymd(data_tmp$`Onset date`) 
      Rcp_case_id <- data_tmp$case_id[which.max(date_tmp)]
      Don_case_id <- data_tmp$case_id[which.min(date_tmp)]
      df$Recipient[i] <<- Rcp_case_id
      df$Donor[i] <<- Don_case_id

      Rcp_lab_id <- df_raw_data %>% filter(`HK case no.` == Rcp_case_id) %>% arrange(desc(`Genome coverage`)) %>% .[1,] %>% .$Sample
      Don_lab_id <- df_raw_data %>% filter(`HK case no.` == Don_case_id) %>% arrange(desc(`Genome coverage`)) %>% .[1,] %>% .$Sample

      read_rcp <- read_delim(paste0("../../2020-09-01_COVID_NGS_pipeline/FastQC/", Rcp_lab_id, "_bam.readcount.mpileup.txt"), "\t")
      read_don <- read_delim(paste0("../../2020-09-01_COVID_NGS_pipeline/FastQC/", Don_lab_id, "_bam.readcount.mpileup.txt"), "\t")
      read_rcp <- read_rcp[!duplicated(read_rcp$pos),]
      read_don <- read_don[!duplicated(read_don$pos),]
      read_rcp <- read_rcp %>% filter(!pos %in% masked_sites)
      read_don <- read_don %>% filter(!pos %in% masked_sites)

      var_list <- mclapply(seq_len(nrow(read_don)), function(ii){
          x <- read_don[ii,]
          depth <- as.numeric(x[["depth"]])
          if(depth==0){
              return(NA)
          }
          count_tmp <- as.numeric(c(x[["acount"]], x[["ccount"]],
                                  x[["gcount"]], x[["tcount"]]))
          tmp <- c("A", "C", "G", "T")
          idx_tmp <- which.max(count_tmp)
          max_nt <- tmp[idx_tmp]
          # check <- count_tmp/depth>=var_calling_threshold & count_tmp/depth <= 0.95 & depth>=100
          check <- (count_tmp/depth >= var_calling_threshold) & (depth>=100)
          check[which.max(count_tmp)] <- FALSE
          if(any(check)){
              # print(ii)
              return(paste0(x$pos, "-", max_nt, "-", tmp[check], "-", x[["depth"]], "-", count_tmp[check]))
          } else {
              return(NA)
          }
      }, mc.cores = mc.cores)

      var_list <- unlist(var_list)
      var_list <- var_list[!is.na(var_list)]
      if(length(var_list)==0){
        return(0)
      }

      # print(length(var_list))
      Position <- sapply(var_list, function(x){
          as.numeric(strsplit(x, "-")[[1]][1])
      })
      Refernece <- sapply(var_list, function(x){
          strsplit(x, "-")[[1]][2]
      })
      Variant <- sapply(var_list, function(x){
          strsplit(x, "-")[[1]][3]
      })
      Depth_donor <- sapply(var_list, function(x){
          as.numeric(strsplit(x, "-")[[1]][4])
      })
      Variant_reads_donor <- sapply(var_list, function(x){
          as.numeric(strsplit(x, "-")[[1]][5])
      })

      data_tmp <- tibble(Pair = id, var_calling_threshold = var_calling_threshold, Position = Position, Refernece = Refernece, Variant = Variant, Depth_donor = Depth_donor, Variant_reads_donor = Variant_reads_donor)

      data_tmp$Depth_recipient <- sapply(seq_len(nrow(data_tmp)), function(ii){
          x <- read_rcp %>% filter(pos == data_tmp$Position[ii])
          as.numeric(x$depth[1])
      })
      data_tmp$Variant_reads_recipient <- sapply(seq_len(nrow(data_tmp)), function(ii){
          x <- read_rcp %>% filter(pos == data_tmp$Position[ii])
          count_tmp <- as.numeric(c(x[["acount"]], x[["ccount"]],
                                  x[["gcount"]], x[["tcount"]]))
          count_tmp[which(c("A", "C", "G", "T") == data_tmp$Variant[ii])]
      })
      
      data_tmp <- data_tmp %>% filter(Depth_recipient>=100)
  
      donor_freqs_observed <- as.data.frame(data_tmp$Variant_reads_donor/data_tmp$Depth_donor)
      n_variants <- nrow(data_tmp)
      recipient_total_reads <- as.data.frame(data_tmp$Depth_recipient) # read.table(args[2])
      recipient_var_reads_observed <- as.data.frame(data_tmp$Variant_reads_recipient) # read.table(args[3])
      recipient_freqs_observed <- recipient_var_reads_observed/recipient_total_reads

      log_likelihood_function <- generate_log_likelihood_function_approx(donor_freqs_observed, recipient_total_reads, recipient_var_reads_observed, Nb_min, Nb_max, var_calling_threshold, confidence_level, n_variants)

      # log_likelihood_function_exact <- generate_log_likelihood_function_exact(donor_freqs_observed, recipient_total_reads, recipient_var_reads_observed, Nb_min, Nb_max, var_calling_threshold, confidence_level, n_variants)

      log_likelihood_function <- restrict_log_likelihood(log_likelihood_function, Nb_min, Nb_max)
      df$bottleneck_size[i] <<- return_bottleneck_size(log_likelihood_function, Nb_min, Nb_max)
      df$CI_index_lower[i] <<- return_CI_lower(log_likelihood_function, Nb_min, Nb_max)
      df$CI_index_upper[i] <<- return_CI_upper(log_likelihood_function, Nb_min, Nb_max)
      log_likelihood_matrix[i,] <<- log_likelihood_function
      df_trans_pair_reads <<- bind_rows(df_trans_pair_reads, data_tmp)
  })

  df$pair <- as.character(df$pair)
  df$Pair <- paste0(df$Donor, " to ", df$Recipient, " (", gsub("Cluster_", "", df$pair), ")")
  df_bottleneck_all <- bind_rows(df_bottleneck_all, df)
}

df_bottleneck_all <- left_join(df_bottleneck_all, df_pair %>% dplyr::select(pair, date_diff) %>% unique())
write_csv(df_bottleneck_all, "../../results/Bottleneck_and_SNVs/bottleneck_1day.csv")


