library(NELSI)
library(tidyverse)
library(lubridate)
library(ggtree)
library(scales)
library(Rmisc)
library(broom)
source("get_clades.R")


trees <- tail(read.nexus("data/wave-1/final-data/local-detection-lag/n238-hk-seqs-alignment-renamed-ref.trees"), 1000)

# Transmission lineages connected to importations are those with the following regular expression tag:
tag <- "Locally|Import_related"

#allocate memory
detection_lag_credible_summary <- data.frame(paste0("Tree ", 1:length(trees))) %>%
  as_tibble() %>%
  transmute(
    tree = paste0..Tree....1.length.trees..,
    tmrca = as_date("2020-01-01"),
    detection_lag = as.numeric(1:length(trees)),
    num_samples = as.numeric(1:length(trees))
  ) %>%
  nest(tmrca, detection_lag, num_samples)

for (tr in 1:length(trees)) {
  hcc_tree <- trees[[tr]]
  youngest_date <- max(decimal_date(ymd(gsub(".+[|]", "", hcc_tree$tip.label))))
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
  tree_summary[, 4] <- num_samples
  tree_summary[, 5] <- sapply(clades, function(x) paste0(x, collapse = ";"))

  #### calculate detection lag per tree
  nested_tree_summary_tbl <- tree_summary %>%
    as_tibble() %>%
    mutate(
      tmrca = as_date(date_decimal(as.numeric(tmrca))),
      first_sample_date = as_date(date_decimal(as.numeric(first_sample_date))),
      last_sample_date = as_date(date_decimal(as.numeric(last_sample_date)))
    ) %>%
    mutate(detection_lag = as.numeric(first_sample_date - tmrca)) %>%
    select(tmrca, detection_lag, num_samples) %>%
    nest(tmrca, detection_lag, num_samples)

  detection_lag_credible_summary[tr, 2] <- nested_tree_summary_tbl[1]

  print(detection_lag_credible_summary[tr, ])
}

detection_lag_credible_summary %>%
  unnest(cols = data) %>%
  mutate(month = month(tmrca)) %>% 
  filter(month == 1) %>%
  pull(detection_lag) %>%
  CI(0.95)

detection_lag_credible_summary %>%
  unnest(cols = data) %>%
  mutate(month = month(tmrca)) %>% 
  group_by(month) %>%
  dplyr::summarise(mean = mean(detection_lag),
                   median = median(detection_lag))

detection_lag_credible_summary %>%
  unnest(cols = data) %>%
  mutate(month = month(tmrca)) %>%
  filter(month == 4) %>% pull(detection_lag) %>% HDInterval::hdi()

  ggplot() +
  geom_histogram(aes(x = detection_lag, y = ..density..), fill = '#dedede', colour = "black", binwidth = 1)
  
  arrange(tmrca)

C <- detection_lag_credible_summary %>%
  unnest(cols = data) %>%
  group_by(tree) %>%
  arrange(tmrca) %>%
  slice(1L) %>%
  ggplot() +
  geom_histogram(aes(x = tmrca, y = ..density..), fill = '#dedede', colour = "black", binwidth = 1) +
  scale_x_date("Date of earliest introduction",
               date_breaks = "14 days", 
               date_labels = "%d %b") + 
  scale_y_continuous("Proportion of posterior density", expand = c(0,0), limits = c(0,0.10), breaks = seq(0,0.1, by = 0.02)) +
  theme_classic() +
  theme(aspect.ratio = 1)



detection_lag_credible_summary %>%
  unnest(cols = data) %>%
  group_by(tree) %>%
  arrange(tmrca) %>%
  slice(1L) %>% 
  mutate(tmrca = decimal_date(tmrca)) %>%
  pull(tmrca) %>% HDInterval::hdi()

date_decimal(2020.000)



 
##Correlation test between detection lag and TMRCA of local lineages
corr_detection <- detection_lag_credible_summary %>%
  unnest(cols = data) %>%
  mutate(tmrca = decimal_date(tmrca)) %>%
  select(tmrca, detection_lag)

cor.test(corr_detection$tmrca, corr_detection$detection_lag,
           method = "spearman") %>%
    tidy()


###Plot sample of detection lag posterior
detection_lag_sample <- detection_lag_credible_summary %>%
  unnest(cols = data) %>%
  mutate(tree_number = parse_number(as.character(tree))) %>%
  sample_n(size = 1000, replace = FALSE)

detection_lag_credible_summary %>%
  unnest(cols = data) %>%
  pull(num_samples) %>% as.numeric() %>% median()




detection_lag_sample


D <- ggplot(detection_lag_sample) +
  geom_jitter(aes(x = tmrca, y = detection_lag), 
              height = 7, 
              width = 7, 
              size = 1, 
              shape = 21,
              alpha = 0.6, 
              fill = "#dedede") +
  geom_smooth(method = lm, aes(x = tmrca, y = detection_lag), 
              se = F, 
              size = 0.7, 
              colour = "black", 
              alpha = 0.3) +
  scale_y_continuous("Detection lag (days)", 
                     expand = c(0, 0), 
                     breaks = seq(0, 70, by = 10), 
                     limits = c(0, 70), 
                     oob = squish) +
  scale_x_date("TMRCA",
    date_breaks = "21 days", 
    date_labels = "%d %b",
    limits = c(date("2020-01-01"), 
               date("2020-04-05"))
  ) +
  theme_classic() +
  theme(aspect.ratio = 1)


row <- cowplot::plot_grid(B, C, D, labels = c('B', 'C', 'D'), nrow = 1, align = 'hv')
