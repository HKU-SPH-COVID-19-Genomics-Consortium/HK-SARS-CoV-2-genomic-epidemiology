library(tidyverse)
library(readxl)

cases_3 <- readLines("../../data/Phylogeography/clade3_headers.txt")
cases_4 <- readLines("../../data/Phylogeography/clade4_headers.txt")
df_rst <- tibble(bind_rows(tibble(seq_name = cases_3, wave = 3), tibble(seq_name = cases_4, wave = 4)))

df_metadata <- read_excel("../../data/metadata_2021-04-30.xlsx", col_types = "text")

data_cluster <- df_metadata %>% select(case.id, Cluster:`Quaternary generation 1`) %>% filter(Cluster == "Y") %>% pivot_longer(`Cluster setting (primary 1)`:`Quaternary generation 1`) %>% filter(!is.na(value))

df_primary_type <- lapply(df_rst$seq_name, function(x){
	case_id <- strsplit(x, "/", fixed = T)[[1]][1]
	case_id <- strsplit(case_id, "_", fixed = T)[[1]][2]
	pri_type <- data_cluster$value[data_cluster$case.id == case_id][1]
	if(grepl("_sec_", pri_type) | grepl("_tert_", pri_type)){
		pri_type <- strsplit(pri_type, "_")[[1]][3]
	} else{
		pri_type <- strsplit(pri_type, "_")[[1]][2]
	}
	return(tibble(seq_name = x, case_id = case_id, type = pri_type))
})
df_primary_type <- bind_rows(df_primary_type)

sort(table(df_primary_type$type))
df_primary_type$type[df_primary_type$type %in% c("friends", "meal", "dating", "MK", "GreenRiver", "BunKeeCongee", "Bldg", "Windsor", "mahjong", "TuenMun", "VictoriaHarbour", "gathering", "hotel", "tutorial", "TST", "table", "Sunfat", "restaurant", "luckydragon", "Delux", "Dating", "boattrip", "boardinghouse")] <- "Social"
df_primary_type$type[df_primary_type$type %in% c("roommate", "fam")] <- "Family/Roommate"
df_primary_type$type[df_primary_type$type %in% c("maid quarter", "work", "TMH", "taxi", "school")] <- "Work"
df_primary_type$type[df_primary_type$type %in% c("RCHE", "RCHD")] <- "RCHE/RCHD"
df_primary_type$type[df_primary_type$type %in% c("GPconsultation")] <- "Unknown"
df_primary_type$type[is.na(df_primary_type$type)] <- "Unknown"

df_primary_type$type <- str_to_title(df_primary_type$type)
sort(table(df_primary_type$type))

# data_cluster$value[grep("GPconsultation", data_cluster$value)]
df_out <- left_join(df_rst, df_primary_type, "seq_name") %>% unique()
write_csv(df_out, "../../results/Phylogeography/transmission_type.csv")
