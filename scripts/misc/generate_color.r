library(tidyverse)
library(RColorBrewer)

data <- read_csv("../../data/seq_global_analysis/metadata_pango_ref_1279.csv")

lineages <- data$`Pango lineage`
mycolors_list <- colorRampPalette(c("#377eb8","#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628"))(length(lineages))
write_csv(tibble(lineage = sort(lineages), color = mycolors_list), "../../data/misc/color_pango_ref.csv")
