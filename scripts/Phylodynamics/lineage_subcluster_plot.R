library(NELSI)
library(tidyverse)
library(lubridate)
library(ggtree)
source("get_clades.R")

setwd("~/Dropbox/GitHub/covid-19-genomics")

#Highest-clade-creibility tree
hcc_tree <- read.nexus('data/wave-1/final-data/summarise-importations/alignment-ref-trimmed-epi-dates.MCC.tree')

hk_lineage_clusters <- read_csv(file = "data/hk_lineage_data.csv")

read_csv(file = "data/wave-1/data/case_data.csv") %>%
  filter(cluster.id != 0) %>% group_by(cluster.id) %>%
  dplyr::mutate(n = n()) %>%
  ungroup() %>%
  group_by(cluster.category) %>%
  dplyr::summarise(median = median(n))
                       
                     
#Some summary data from the tree
youngest_date <- max(decimal_date(ymd(gsub('.+[|]', '', hcc_tree$tip.label))))
tips <- hcc_tree$tip.label
node_times <- youngest_date - allnode.times(hcc_tree, reverse = F)

#Transmission lineages and importations are those with the following regular expression tag:
tag <- 'Locally|Import_related'

# Find monophyletic groups and singletons with the tag above
clades <- find.monophyletic(hcc_tree, tag, include.singletons = T)

# Get closest sister taxon, if the sister is a clade take the oldest tip  
first_sisters <- vector()
first_sister_dates <- vector()
for(i in 1:length(clades)){
  sister_temp <- hcc_tree$tip.label[find.sister(hcc_tree, which(hcc_tree$tip.label %in% clades[[i]]))]
  sister_dates <- decimal_date(ymd(gsub('.+[|]', '', sister_temp)))
  first_sisters[i] <- sister_temp[which.min(sister_dates)]
  first_sister_dates[i] <- min(sister_dates)
}

# Get tmrcas of clades
clades_indexed <- lapply(clades, function(x) which(hcc_tree$tip.label %in% x))
tmrcas <- sapply(clades_indexed,
                 function(x)
                   node_times[get.mrca(hcc_tree, x)])

# Find oldest and youngest sampling times for each clade
first_samples <- sapply(clades_indexed, function(x) min(node_times[x]))
last_samples <- sapply(clades_indexed, function(x) max(node_times[x]))

# Number of samples per clade
num_samples <- sapply(clades, function(x) length(x))

#Make summary matrix
tree_summary <- matrix(NA, length(clades), 7)
colnames(tree_summary) <- c('tmrca', 'first_sample_date',
                            'last_sample_date', 'num_samples',
                            'samples', 'first_sister_tip',
                            'first_sister_date')

tree_summary[, 1] <- tmrcas
tree_summary[, 2] <- first_samples
tree_summary[, 3] <- last_samples
tree_summary[, 4] <- num_samples
tree_summary[, 5] <- sapply(clades, function(x) paste0(x, collapse = ';'))
tree_summary[, 6] <- first_sisters
tree_summary[, 7] <- first_sister_dates


tree_summary_clades <- tree_summary %>%
  as_tibble() %>%
  mutate(
    tmrca = as_date(date_decimal(as.numeric(tmrca))),
    first_sample_date = as_date(date_decimal(as.numeric(first_sample_date))),
    last_sample_date = as_date(date_decimal(as.numeric(last_sample_date)))
  ) %>%
  arrange(tmrca) %>%
  mutate(
    clade = 1:nrow(.),
    case = str_extract_all(samples, "HK\\d+"),
    sister_case = str_extract(first_sister_tip, "HK\\d+")
  ) %>%
  unnest(case) %>%
  mutate(
    hk_case_no = parse_number(case),
    sister_case = parse_number(sister_case) 
  ) %>%
  distinct() %>%
  select(hk_case_no, sister_case, clade) %>%
  unite(case_no, hk_case_no, sister_case) %>%
  mutate(case_no = str_extract_all(case_no, "\\d+")) %>%
  unnest(case_no) %>%
  dplyr::mutate(case_no = as.numeric(case_no)) %>%
  distinct()
    
  
ct_cluster_tbl <- read_csv(file = "data/wave-1/data/case_data.csv") %>% janitor::clean_names() %>%
  dplyr::mutate(cluster_id = case_when(case_no == 227 ~ 0, #correct contact tracing misclassificaiton
                                       TRUE ~ cluster_id),
                case_new_classification = case_when(case_classification == "Imported" ~ "Travel overseas",
                                                    case_classification == "Local case" ~ "Local unknown source",
                                                    TRUE ~ "Contact with a confirmed case"))


binded_case_data <- as_tibble()
tree_summary_clades_list <- tree_summary_clades %>% select(clade) %>% distinct()

for (i in seq_along(tree_summary_clades_list$clade)) {
  
  cluster_link <- ct_cluster_tbl %>%
    left_join(tree_summary_clades, by = "case_no") %>% filter(clade == tree_summary_clades_list$clade[i]) %>%
    dplyr::select(cluster_id, clade) %>% distinct() %>% filter(cluster_id != 0)
  
  new_case_data <- ct_cluster_tbl %>%
    left_join(tree_summary_clades, by = "case_no") %>% filter(clade == tree_summary_clades_list$clade[i] | cluster_id %in% cluster_link$cluster_id) %>%
    dplyr::mutate(metacluster_id = tree_summary_clades_list$clade[i], 
                  seqs = case_when(!is.na(clade) ~ "Y",
                                   TRUE ~ "N")) %>% 
    select(-clade) %>% distinct()
  
  binded_case_data <- bind_rows(binded_case_data, new_case_data)
  
}

#####Further account for PANGO lineages within clade designations
hk_lineage <- hk_lineage_clusters %>% transmute(case_no = hk_case_no, lineage = lineage)
binded_lineage_case_data <- as_tibble()
metacluster_lineage <- binded_case_data %>% left_join(hk_lineage, by = "case_no") %>% group_by(metacluster_id, lineage) %>% dplyr::count() %>% filter(!is.na(lineage))

for (i in seq_along(metacluster_lineage$metacluster_id)) {
  
  seqs_lineage_case_data <- binded_case_data %>%  left_join(hk_lineage, by = "case_no") %>% filter(metacluster_id == metacluster_lineage$metacluster_id[i], lineage == metacluster_lineage$lineage[i]) %>% dplyr::mutate(lineage_metacluster_id = i) 
  seqs_lineage_case_data_link <- seqs_lineage_case_data %>% dplyr::select(cluster_id, lineage) %>% distinct() %>% filter(cluster_id != 0)
  unseq_lineage_case_data <- binded_case_data %>% left_join(hk_lineage, by = "case_no") %>% 
    filter(metacluster_id == metacluster_lineage$metacluster_id[i], cluster_id %in% seqs_lineage_case_data_link$cluster_id, seqs == "N") %>% dplyr::mutate(lineage_metacluster_id = i)
  combined_lineage_case_data <- bind_rows(seqs_lineage_case_data, unseq_lineage_case_data)
  binded_lineage_case_data <- bind_rows(binded_lineage_case_data, combined_lineage_case_data)
  
}

#renumber lineages excluding singles
binded_lineage_case_data_factor <- binded_lineage_case_data %>%
  distinct() %>%
  group_by(lineage_metacluster_id) %>%
  dplyr::mutate(n = n()) %>%
  ungroup() %>%
  filter(n > 1) %>%
  dplyr::mutate(id = as_factor(group_indices(. , lineage_metacluster_id))) %>%
  dplyr::mutate(onset_date = dmy(onset_date),
         epi_date = dmy(epi_date),
         confirm_date = dmy(confirm_date))

cluster_count_size <- binded_lineage_case_data_factor %>%
  group_by(id) %>%
  slice(1L) %>%
  pull(n)

cluster_count_size_position <- binded_lineage_case_data_factor %>%
  group_by(id) %>%
  arrange(desc(confirm_date)) %>%
  slice(1L) %>%
  dplyr::mutate(confirm_date = confirm_date+10) %>%
  pull(confirm_date)



##PLOT OF MCC TREE LINEAGE AND NELSI CLUSTER
ggplot(binded_lineage_case_data_factor) +
  geom_jitter(aes(x = confirm_date, y = id, fill = lineage), shape = 21, height = 0.1, width = 2, size = 1.5, alpha = 0.9) +
  annotate("text", x = as.Date(cluster_count_size_position), y = 1:length(cluster_count_size), label = cluster_count_size, size = 2.5, alpha = 0.8) +
  scale_x_date("Report date",
               date_breaks = "21 days", 
               date_labels = "%d %b", 
               limits = as.Date(c("2020-01-19","2020-04-27"))) +
  theme_classic() +
  theme(legend.position = 'right',
        panel.grid.major.y = element_line(size = 0.1)) +
  scale_fill_viridis_d() +
  labs(fill = "Lineage",
      y = "Cluster")


ggsave(filename = "output/final-figures/singles/lineages.pdf", height = 10, width = 4)

binded_lineage_case_data_factor %>%
  filter(id == "5") %>% view()
  group_by(cluster_id) %>% dplyr::count()

############ Calculate  lineage from the posterior set of 1000 trees

trees <- tail(read.nexus("data/wave-1/final-data/summarise-importations/alignment-ref-trimmed-epi-dates.trees"), 1000)

ct_cluster_tbl <- read_csv(file = "data/wave-1/data/case_data.csv") %>% janitor::clean_names() %>%
  dplyr::mutate(cluster_id = case_when(case_no == 227 ~ 0, 
                                       TRUE ~ cluster_id),
                case_new_classification = case_when(case_classification == "Imported" ~ "Imported case",
                                                    TRUE ~ "Local case"))

tag <- 'Locally|Import_related'


#allocate memory
lineage_count_credible_summary <- data.frame(paste0("Tree ", 1:length(trees))) %>%
  as_tibble() %>%
  transmute(
    tree = paste0..Tree....1.length.trees..,
    n = as.numeric(1:length(trees))
  ) %>%
  nest(n)

##Loop to calculate lineage size distirbution across 1000 posteior annotated trees
for (tr in 1:length(trees)) {
  hcc_tree <- trees[[tr]]
  youngest_date <- max(decimal_date(ymd(gsub('.+[|]', '', hcc_tree$tip.label))))
  tips <- hcc_tree$tip.label
  node_times <- youngest_date - allnode.times(hcc_tree, reverse = F)
  clades <- find.monophyletic(hcc_tree, tag, include.singletons = T)
  
  first_sisters <- vector()
  first_sister_dates <- vector()
  for(i in 1:length(clades)){
    sister_temp <- hcc_tree$tip.label[find.sister(hcc_tree, which(hcc_tree$tip.label %in% clades[[i]]))]
    sister_dates <- decimal_date(ymd(gsub('.+[|]', '', sister_temp)))
    first_sisters[i] <- sister_temp[which.min(sister_dates)]
    first_sister_dates[i] <- min(sister_dates)
  }
  
  clades_indexed <- lapply(clades, function(x) which(hcc_tree$tip.label %in% x))
  tmrcas <- sapply(clades_indexed,
                   function(x)
                     node_times[get.mrca(hcc_tree, x)])
  
  # Find oldest and youngest sampling times for each clade
  first_samples <- sapply(clades_indexed, function(x) min(node_times[x]))
  last_samples <- sapply(clades_indexed, function(x) max(node_times[x]))
  
  # Number of samples per clade
  num_samples <- sapply(clades, function(x) length(x))
  
  #Make summary matrix
  tree_summary <- matrix(NA, length(clades), 7)
  colnames(tree_summary) <- c('tmrca', 'first_sample_date',
                              'last_sample_date', 'num_samples',
                              'samples', 'first_sister_tip',
                              'first_sister_date')
  
  tree_summary[, 1] <- tmrcas
  tree_summary[, 2] <- first_samples
  tree_summary[, 3] <- last_samples
  tree_summary[, 4] <- num_samples
  tree_summary[, 5] <- sapply(clades, function(x) paste0(x, collapse = ';'))
  tree_summary[, 6] <- first_sisters
  tree_summary[, 7] <- first_sister_dates
  
  
  tree_summary_clades <- tree_summary %>%
    as_tibble() %>%
    mutate(
      tmrca = as_date(date_decimal(as.numeric(tmrca))),
      first_sample_date = as_date(date_decimal(as.numeric(first_sample_date))),
      last_sample_date = as_date(date_decimal(as.numeric(last_sample_date)))
    ) %>%
    arrange(tmrca) %>%
    mutate(
      clade = 1:nrow(.),
      case = str_extract_all(samples, "HK\\d+"),
      sister_case = str_extract(first_sister_tip, "HK\\d+")
    ) %>%
    unnest(case) %>%
    mutate(
      hk_case_no = parse_number(case),
      sister_case = parse_number(sister_case) 
    ) %>%
    distinct() %>%
    select(hk_case_no, sister_case, clade) %>%
    unite(case_no, hk_case_no, sister_case) %>%
    mutate(case_no = str_extract_all(case_no, "\\d+")) %>%
    unnest(case_no) %>%
    dplyr::mutate(case_no = as.numeric(case_no)) %>%
    distinct()
  
  
  binded_case_data <- as_tibble()
  tree_summary_clades_list <- tree_summary_clades %>% select(clade) %>% distinct()
  
  for (j in seq_along(tree_summary_clades_list$clade)) {
    
    cluster_link <- ct_cluster_tbl %>%
      left_join(tree_summary_clades, by = "case_no") %>% filter(clade == tree_summary_clades_list$clade[j]) %>%
      dplyr::select(cluster_id, clade) %>% distinct() %>% filter(cluster_id != 0)
    
    new_case_data <- ct_cluster_tbl %>%
      left_join(tree_summary_clades, by = "case_no") %>% filter(clade == tree_summary_clades_list$clade[j] | cluster_id %in% cluster_link$cluster_id) %>%
      dplyr::mutate(metacluster_id = tree_summary_clades_list$clade[j], 
                    seqs = case_when(!is.na(clade) ~ "Y",
                                     TRUE ~ "N")) %>% 
      select(-clade) %>% distinct()
    
    binded_case_data <- bind_rows(binded_case_data, new_case_data)
    
  }
  
  hk_lineage <- hk_lineage_clusters %>% transmute(case_no = hk_case_no, lineage = lineage)
  binded_lineage_case_data <- as_tibble()
  metacluster_lineage <- binded_case_data %>% left_join(hk_lineage, by = "case_no") %>% group_by(metacluster_id, lineage) %>% dplyr::count() %>% filter(!is.na(lineage))
  
  for (k in seq_along(metacluster_lineage$metacluster_id)) {
    
    seqs_lineage_case_data <- binded_case_data %>%  left_join(hk_lineage, by = "case_no") %>% filter(metacluster_id == metacluster_lineage$metacluster_id[k], lineage == metacluster_lineage$lineage[k]) %>% dplyr::mutate(lineage_metacluster_id = k) 
    seqs_lineage_case_data_link <- seqs_lineage_case_data %>% dplyr::select(cluster_id, lineage) %>% distinct() %>% filter(cluster_id != 0)
    unseq_lineage_case_data <- binded_case_data %>% left_join(hk_lineage, by = "case_no") %>% 
      filter(metacluster_id == metacluster_lineage$metacluster_id[k], cluster_id %in% seqs_lineage_case_data_link$cluster_id, seqs == "N") %>% dplyr::mutate(lineage_metacluster_id = k)
    combined_lineage_case_data <- bind_rows(seqs_lineage_case_data, unseq_lineage_case_data)
    binded_lineage_case_data <- bind_rows(binded_lineage_case_data, combined_lineage_case_data)
    
  }
  
  binded_lineage_case_data_factor <- binded_lineage_case_data %>%
    distinct() %>%
    group_by(lineage_metacluster_id) %>%
    dplyr::mutate(n = n()) %>%
    ungroup() %>%
    filter(n > 1) %>%
    dplyr::mutate(id = as_factor(group_indices(. , lineage_metacluster_id))) %>%
    dplyr::mutate(onset_date = dmy(onset_date),
                  epi_date = dmy(epi_date),
                  confirm_date = dmy(confirm_date))
  
  
  nested_binded_lineage_case_data_factor <- binded_lineage_case_data_factor %>%
    group_by(id) %>%
    slice(1L) %>%
    ungroup() %>%
    dplyr::select(n) %>%
    nest(data = c(n))


lineage_count_credible_summary[tr, 2] <- nested_binded_lineage_case_data_factor[1]

print(lineage_count_credible_summary[tr, ])

}
  

lineage_count_credible_summary %>%
  unnest(cols = data) %>% pull(n) %>% range()
  
  
  arrange(desc(n)) %>% filter(n != 107) %>% pull(n) %>% HDInterval::hdi()


#group_by(tree) %>% dplyr::count() %>% pull(n) %>% max()

pull(n) %>% 
  
  
  
  
  ggplot() +
  geom_histogram(aes(x = n, y = ..density..), fill = '#dedede', colour = "black", binwidth = 0.25) +
  scale_x_log10()
  scale_x_date("Date of earliest introduction",
               date_breaks = "14 days", 
               date_labels = "%d %b") + 
  scale_y_continuous("Proportion of posterior density", expand = c(0,0), limits = c(0,0.10), breaks = seq(0,0.1, by = 0.02)) +
  theme_classic() +
  theme(aspect.ratio = 1)