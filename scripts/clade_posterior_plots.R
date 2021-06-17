library(tidyverse)
library(lubridate)
library(scales)
set.seed(1)


##### Due to the large size of the results file, it is not available via GitHub but can be downloaded from the link below:  https://www.dropbox.com/s/es7lhheemlpgm2d/credible_tree_summary_sample_8000.csv?dl=0

credible_tree_summary <- read_csv("results/credible_tree_summary_sample_8000.csv")

credible_tree_summary %>%
  mutate(wave = case_when(first_sample_date <= "2020-02-22" ~ "Wave 1",
                          first_sample_date >= "2020-02-23" & first_sample_date <= "2020-05-12" ~ "Wave 2",
                          first_sample_date >= "2020-05-13" & first_sample_date <= "2020-09-29" ~ "Wave 3",
                          first_sample_date >= "2020-09-29" & first_sample_date <= "2021-01-26" ~ "Wave 4",
                          TRUE ~ NA_character_),
         clade = case_when(num_travel_cases == 0 & num_community_cases == 1 ~ "Single local",
                           num_travel_cases > 0 & num_community_cases == 1 ~ "Singleton",
                           num_travel_cases > 0 & num_community_cases > 1 ~ "Community lineage",
                           num_travel_cases == 0 & num_community_cases != 0 ~ "Community lineage",
                           num_travel_cases > 0 & num_community_cases == 0 ~ "Single import",
                           TRUE ~ NA_character_)) %>% 
  filter(wave == "Wave 3" | wave == "Wave 4",
         clade == "Community lineage") %>%
  pull(effective_delay) %>%
  HDInterval::hdi()


### Number of lineages per wave
credible_tree_summary %>%
  filter(tree == "Tree 1") %>%
  mutate(wave = case_when(first_sample_date <= "2020-02-22" ~ "Wave 1",
                          first_sample_date >= "2020-02-23" & first_sample_date <= "2020-05-12" ~ "Wave 2",
                          first_sample_date >= "2020-05-13" & first_sample_date <= "2020-09-29" ~ "Wave 3",
                          first_sample_date >= "2020-09-29" & first_sample_date <= "2021-01-26" ~ "Wave 4",
                          TRUE ~ NA_character_),
         singleton = case_when(num_samples == 1 & num_community_cases != 0 ~ "Singleton",
                               num_travel_cases > 0 & num_community_cases == 0 ~ "Singleton",
                               TRUE ~ "Lineage")) %>%
  group_by(wave, singleton) %>%
  dplyr::count() %>%
  spread(key = singleton, value = n) %>%
  mutate(total = Lineage + Singleton)


## TMRCA of the five earliest lineage introductions
lineage_n19 <- credible_tree_summary %>%
  filter(str_detect(samples, "HKcase_24\\|")) %>%
  mutate(lineage = "n=19")

lineage_n7 <- credible_tree_summary %>%
  filter(str_detect(samples, "HKcase_42\\|")) %>%
  mutate(lineage = "n=7")

lineage_n92 <- credible_tree_summary %>%
  filter(str_detect(samples, "HKcase_192\\|")) %>%
  mutate(lineage = "n=92")

lineage_n29 <- credible_tree_summary %>%
  filter(str_detect(samples, "HKcase_923\\|")) %>%
  mutate(lineage = "n=29")

lineage_n6 <- credible_tree_summary %>%
  filter(str_detect(samples, "HKcase_229\\|")) %>%
  mutate(lineage = "n=6")



A <- bind_rows(lineage_n19, lineage_n7, lineage_n92, lineage_n29, lineage_n6) %>%
  ggplot() +
  geom_density(aes(x=tmrca, colour = lineage), adjust = 20, size = 0.8) +
  scale_x_date("Lineage TMRCA",
               date_breaks = "1 month", 
               date_labels = "%b/%y",
               expand = c(0,0),
               limits = c(date("2019-12-15"), date("2020-05-12"))) + 
  scale_y_continuous("Model density",
                     limits = c(0,0.15),
                     breaks = seq(0,0.15, by = 0.03),
                     expand = c(0,0)) +
  theme_classic() +
  theme(aspect.ratio = 1, 
        legend.position = c(0.85, 0.85), 
        legend.title = element_blank()) +
  ggsci::scale_color_jco()


#### Number of samples per lineage 
B <- credible_tree_summary %>%
  filter(tree == "Tree 1") %>%
  mutate(num_samples = case_when(num_samples == 100 ~ 902,
                                 num_samples == 65 ~ 552,
                                 TRUE ~ num_samples),
         mean_date = as_date((as.numeric(first_sample_date) + as.numeric(last_sample_date))/2)) %>%
  ggplot() +
  geom_jitter(aes(y = num_samples, 
                  x = first_sample_date, 
                  size = num_samples),
              alpha = 0.8,
              shape = 21, 
              fill = "#dedede",
              height = 0.01) +
  scale_y_continuous("Samples per lineage",
                expand = c(0.01,0), 
                limits = c(NA,1000),
                trans = "log10",
                oob = squish) +
  scale_x_date("Date of lineage detection",
               date_breaks = "3 month", 
               date_labels = "%b/%y",
               expand = c(0.1,0.1)) +
  theme_classic() +
  theme(aspect.ratio = 1,
        legend.position = 'none') +
  annotation_logticks()


##Local Lineage duration vs delay in detection

credible_tree_summary_sample <- credible_tree_summary %>% 
  filter(num_samples > 1,
         num_travel_cases == 0) %>%
  arrange(desc(num_samples)) %>%
  mutate(num_samples = case_when(num_samples == 100 ~ 902,
                                 num_samples == 65 ~ 552,
                                 TRUE ~ num_samples)) %>%
  sample_n(size = 1000, replace = FALSE)



C <- credible_tree_summary_sample %>%
  mutate(wave = case_when(first_sample_date <= "2020-05-12" ~ "Wave 1+2",
                        first_sample_date >= "2020-05-13" & first_sample_date <= "2020-09-29" ~ "Wave 3",
                        first_sample_date >= "2020-09-29" & first_sample_date <= "2021-01-26" ~ "Wave 4",
                        TRUE ~ NA_character_)) %>%
  ggplot() +
  geom_jitter(aes(x = delay, y = undetected_duration, size = num_samples, fill = wave), 
              height = 5, 
              width = 5, 
              shape = 21,
              alpha = 0.7) +
  geom_smooth(method = lm, aes(x = delay, y = undetected_duration),
              se = T, 
              size = 0.7, 
              colour = "black", 
              alpha = 0.3) +
  scale_x_continuous("Detection lag (days)",
                     limits = c(0,60),
                     breaks = seq(0,60, by = 10),
                     expand = c(0, 0), 
                     oob = squish) +
  scale_y_continuous("Lineage duration (days)",
                     limits = c(0,180),
                     breaks = seq(0,180, by = 30),
                     expand = c(0, 0), 
                     oob = squish) +
  theme_classic() +
  theme(aspect.ratio = 1, 
        legend.position = c(0.85,0.85)) +
  labs(size = "Samples per lineage") +
  scale_fill_manual(values = c('#52cc83', '#5283cc', '#b452cc'))

C

detection_duration_correlation <- credible_tree_summary %>% 
  filter(num_samples > 1,
         num_travel_cases == 0) %>%
  arrange(desc(num_samples)) %>%
  mutate(num_samples = case_when(num_samples == 100 ~ 902,
                                 num_samples == 65 ~ 552,
                                 TRUE ~ num_samples))


cor.test(detection_duration_correlation$effective_delay, detection_duration_correlation$undetected_duration,
         method = "spearman",
         exact = F)

##Plot top row
top_row <- cowplot::plot_grid(A, B, C, labels = c('A', 'B', 'C'), align = 'hv', ncol = 1)


### Detection over time by wave
D <- credible_tree_summary %>%
  mutate(local = str_detect(samples, "local"),
         wave = case_when(first_sample_date >= "2020-01-23" & last_sample_date <= "2020-05-12" ~ "Wave 1+2",
                          first_sample_date >= "2020-05-13" & last_sample_date <= "2020-09-29" ~ "Wave 3",
                          first_sample_date >= "2020-09-30" & last_sample_date <= "2021-01-26" ~ "Wave 4",
                          TRUE ~ NA_character_)) %>%
  filter(local == T,
         effective_delay != 0) %>%
  filter(!is.na(wave)) %>%
  group_by(wave) %>%
  sample_n(size = 1000, replace = F) %>%
  ggplot() +
  geom_jitter(aes(x = as_date(tmrca), y = effective_delay, fill = wave),
              height = 10, 
              width = 10,
              shape = 21,
              alpha = 0.6) +
  geom_smooth(method = lm, aes(x = as_date(tmrca), y = effective_delay, fill = wave, colour = wave), 
              se = T, 
              alpha = 0.5) +
  facet_wrap(~wave, scales = 'free_x', ncol = 1) +
  scale_y_continuous("Detection lag (days)", 
                     limits = c(0,100),
                     breaks = seq(0,100, by = 10),
                     expand = c(0, 0), 
                     oob = squish) +
  scale_x_date("TMRCA",
               date_breaks = "1 month", 
               date_labels = "%b/%y") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = c(0.85,0.85),
        panel.spacing = unit(1.5, "lines"),
        aspect.ratio = 1) +
  scale_color_manual(values = c('#52cc83', '#5283cc', '#b452cc')) +
  scale_fill_manual(values = c('#52cc83', '#5283cc', '#b452cc')) +
  labs(fill = "Wave", colour = "Wave")

