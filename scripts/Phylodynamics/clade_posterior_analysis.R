library(NELSI)
library(tidyverse)
library(lubridate)
library(ggtree)
library(scales)
library(Rmisc)
library(broom)
source("get_clades.R")


##### Due to memory limitations in Rstudio, this script was run in the terminal application on macOS after changeing the limit with ```ulimit -s unlimited```

##### Due to the large size of the posterior tree file, it is not available via GitHub but can be downloaded from the link below:
## https://www.dropbox.com/s/sy4uk3qqe241o4l/global_sub_v2.afa.strict.no.n8000.sample.trees?dl=0

print("Reading Tree File")

trees <- read.nexus("data/global/global_sub_v2.afa.strict.no.n8000.sample.trees")

print("Read Complete")

# Transmission lineages connected to importations are those with the following regular expression tag:
tag <- "HKcase"

#allocate memory
credible_tree_summary <- data.frame(paste0("Tree ", 1:length(trees))) %>%
  as_tibble() %>%
  transmute(
    tree = paste0..Tree....1.length.trees..,
    tmrca = as_datetime("2020-01-01"),
    first_sample_date = as_datetime("2020-01-01"),
    last_sample_date = as_datetime("2020-01-01"),
    num_samples = as.numeric(1:length(trees)),
    samples = as.character(1:length(trees)),
    index = as.numeric(1:length(trees)),
    num_travel_cases = as.numeric(1:length(trees)),
    num_community_cases = as.numeric(1:length(trees)),
    delay = as.numeric(1:length(trees)),
    effective_delay = as.numeric(1:length(trees)),
    known_duration = as.numeric(1:length(trees)),
    undected_duration = as.numeric(1:length(trees))
  ) %>%
  nest(data = c(tmrca, first_sample_date, last_sample_date, num_samples, samples, index, 
                num_travel_cases, num_community_cases, delay, effective_delay, 
                known_duration, undected_duration))


for (tr in 1:length(trees)) {
  hcc_tree <- trees[[tr]]
  
  
  youngest_date <- str_extract_all(hcc_tree$tip.label, "\\d+\\-\\d+\\-\\d+") %>%
    ymd() %>% 
    as_tibble() %>%
    filter(!is.na(value)) %>%
    mutate(value = decimal_date(value)) %>%
    pull(value) %>% 
    max()
  
  tips <- hcc_tree$tip.label
  node_times <- youngest_date - allnode.times(hcc_tree, reverse = F)
  clades <- find.monophyletic(hcc_tree, tag, include.singletons = T)
  
  
  clades_indexed <- lapply(clades, function(x) which(hcc_tree$tip.label %in% x))
  tmrcas <- sapply(
    clades_indexed,
    function(x) {
      node_times[get.mrca(hcc_tree, x)]
    }
  )
  
  first_samples <- sapply(clades_indexed, function(x) min(node_times[x]))
  last_samples <- sapply(clades_indexed, function(x) max(node_times[x]))
  num_samples <- sapply(clades, function(x) length(x))
  
  tree_summary <- matrix(NA, length(clades), 5)
  
  colnames(tree_summary) <- c(
    "tmrca", "first_sample_date",
    "last_sample_date", "num_samples",
    "samples"
  )
  
  tree_summary[, 1] <- tmrcas
  tree_summary[, 2] <- first_samples
  tree_summary[, 3] <- last_samples
  tree_summary[, 4] <- as.numeric(num_samples)
  tree_summary[, 5] <- sapply(clades, function(x) paste0(x, collapse = ";"))
  
  #### calculate detection lag per tree
  # nested_tree_summary_tbl <- 
  
  tree_summary_tbl <- as_tibble(tree_summary) %>% 
    mutate(index = 1:nrow(.),
           tmrca = as_datetime(date_decimal(as.numeric(tmrca))),
           first_sample_date = as_datetime(date_decimal(as.numeric(first_sample_date))),
           last_sample_date = as_datetime(date_decimal(as.numeric(last_sample_date))),
           num_travel_cases = 0,
           num_community_cases = 0,
           num_samples = as.numeric(num_samples))
  
  for (i in seq_along(tree_summary_tbl$index)) {
    
    tree_summary_tbl[i,] <- tree_summary_tbl %>%
      filter(index == i) %>%
      mutate(num_travel_cases = as.numeric(length(unlist(str_extract_all(samples, "\\|Imported\\|")))),
             num_community_cases = num_samples - num_travel_cases)
    
  }
  
  nested_tree_summary_tbl <- tree_summary_tbl %>%
    mutate(delay = as.numeric(first_sample_date) - as.numeric(tmrca),
           effective_delay = case_when(num_travel_cases != 0 ~ 0,
                                       TRUE ~ delay),
           known_duration = case_when(num_community_cases >= 2 ~ as.numeric(last_sample_date) - as.numeric(first_sample_date),
                                      TRUE ~ 0),
           undetected_duration = case_when(num_community_cases >= 2 ~ as.numeric(last_sample_date) - as.numeric(tmrca),
                                           TRUE ~ 0)) %>%
    nest(data = everything())
  
  
  
  credible_tree_summary[tr, 2] <- nested_tree_summary_tbl[1]
  
  print(credible_tree_summary[tr, ])
}

credible_tree_summary %>% 
  unnest(cols = c(data)) %>%
  write_csv(file = "credible_tree_summary_sample_8000.csv")
