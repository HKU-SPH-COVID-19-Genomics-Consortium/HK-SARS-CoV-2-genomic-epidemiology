library(tidyverse)
library(readxl)
library(parallel)

cases_3 <- readLines("../../data/Misc/clade3_headers.txt")
cases_4 <- readLines("../../data/Misc/clade4_headers.txt")
df_rst <- tibble(bind_rows(tibble(seq_name = cases_3, wave = 3), tibble(seq_name = cases_4, wave = 4)))
df_rst$case.id <- sapply(df_rst$seq_name, function(x){
	case_id <- strsplit(x, "/", fixed = T)[[1]][1]
	case_id <- strsplit(case_id, "_", fixed = T)[[1]][2]
})

df_metadata <- read_excel("../../data/metadata_2021-04-30.xlsx", col_types = "text")
table(df_metadata$reclassification)
# df_metadata <- df_metadata %>% filter(grepl("ocal", reclassification)) 
df_metadata <- df_metadata %>% filter(reclassification != "Imported") 
table(df_metadata$reclassification)

df_metadata$wave <- mclapply(df_metadata$case.id, function(x){
	print(x)
	date_t <- df_metadata %>% filter(case.id==x) %>% .$`Report date`
	date_t <- lubridate::ymd(date_t)
	if(length(date_t)==0){
		return(NA)
	}
	if(is.na(date_t)){
		return(NA)
	}
	if(date_t <= lubridate::ymd("2020-02-22")){
		return("Wave_1")
	} else if(date_t <= lubridate::ymd("2020-05-12")){
		return("Wave_2")
	} else if(date_t <= lubridate::ymd("2020-09-29")){
		return("Wave_3")
	} else if(date_t <= lubridate::ymd("2021-01-26")){
		return("Wave_4")		
	}
}, mc.cores = 24)
df_metadata$wave <- unlist(df_metadata$wave)

data_cluster <- df_metadata %>% select(case.id, reclassification, Cluster:`Quaternary generation 1`) %>% filter(Cluster == "Y") %>% pivot_longer(`Cluster setting (primary 1)`:`Quaternary generation 1`) %>% filter(!is.na(value))

data_cluster %>% filter(grepl("care", value))

df_primary_type <- lapply(df_metadata$case.id, function(x){
	pri_type <- data_cluster$value[data_cluster$case.id == x][1]
	if(grepl("_sec_", pri_type) | grepl("_tert_", pri_type) | grepl("_quad_", pri_type) | grepl("_penta_", pri_type)){
		pri_type <- strsplit(pri_type, "_")[[1]][3]
	} else{
		pri_type <- strsplit(pri_type, "_")[[1]][2]
	}
	return(tibble(case.id = x, type = pri_type))
})
df_primary_type <- bind_rows(df_primary_type)

sort(table(df_primary_type$type)) 
## social
df_primary_type$type[df_primary_type$type %in% c("friends", "meal", "dating", "MK", "GreenRiver", "BunKeeCongee", "Bldg", "Windsor", "mahjong", "TuenMun", "VictoriaHarbour", "gathering", "hotel", "tutorial", "TST", "table", "Sunfat", "restaurant", "luckydragon", "Delux", "Dating", "boattrip", "boardinghouse", "港九小輪", "royal garden", "China Secret", "singing", "seaview resort", "Gateway", "danceclub", "Bar and Band", "worship", "karaoke", "Otto", "GlowSpa", "土瓜灣喜宴酒樓", "privatetutor", "gym", "DimSumSquare", "church", "wedding", "friend", "sheltered", "Kagesha2", "Kagesha1", "BirthdayParty", "bar", "松山", "travel", "Superstar", "outback TST", "M&S", "hugging", "gym45", "facial", "ChinaBar", "chatting", "beautyparlour", "Bus", "Coffee", "Travel")] <- "Social"
## Family/Roommate
df_primary_type$type[df_primary_type$type %in% c("roommate", "fam", "麗晶Block6", "Block 8 Kwai Shing West", "Holly Mansion", "星月居", "estate", "care", "星月樓", "泉章居", "LukChuen", "明麗樓", "building", "仁石樓Unit 09", "Roommate", "neighbour", "MeiFoo", "maid", "麗港城第五座", "貴東樓", "Partner")] <- "Family/Roommate"
## Work
df_primary_type$type[df_primary_type$type %in% c("maid quarter", "work", "TMH", "taxi", "school", "Kerry Logistics", "KinWing")] <- "Work"
## RCHE/RCHD
df_primary_type$type[df_primary_type$type %in% c("RCHE", "RCHD")] <- "RCHE/RCHD"
## Unknown
df_primary_type$type[df_primary_type$type %in% c("GPconsultation", "CHP", "contact", "ambulance")] <- "Unknown/Sporadic"
df_primary_type$type[is.na(df_primary_type$type)] <- "Unknown/Sporadic"

df_primary_type$type <- str_to_title(df_primary_type$type)
sort(table(df_primary_type$type))

# data_cluster %>% filter(grepl("contact", value)) %>% .$value
# data_cluster %>% filter(grepl("sheltered", value))

df_seqed <- left_join(df_rst, df_primary_type, "case.id") %>% unique()
write_csv(df_seqed, "../../results/Misc/transmission_type_Sequenced.csv")
df_seqed_stat <- df_seqed %>% group_by(wave) %>% summarise(n_social = sum(type == "Social", na.rm = T), n_non_social = sum(type %in% c("Family/Roommate", "Work", "Rche/Rchd", "Nosocomial", na.rm = T)))

# summary count for all cases
df_all <- df_primary_type %>% left_join(df_metadata %>% select(case.id, wave))
write_csv(df_all, "../../results/Misc/transmission_type_All.csv")
df_all %>% filter(wave %in% c("Wave_3", "Wave_4")) %>% group_by(type, wave) %>% summarise(n = n())  
df_all %>% group_by(wave) %>% summarise(n = n())


df_all_stat <- df_all %>% group_by(wave) %>% summarise(n_social = sum(type == "Social", na.rm = T), n_non_social = sum(type %in% c("Family/Roommate", "Work", "Rche/Rchd", "Nosocomial", na.rm = T)))

# chisq test
chisq.test(df_seqed_stat[1:2, 2:3])
chisq.test(df_all_stat[3:4, 2:3])

