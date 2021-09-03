library("tidyverse")
library("readxl")

# Were the samples of epi-linked people (ext data table 8) sequenced during the same run (to avoid between run variability)? If not, this should be detailed, or a note should be made on why this does/does not affect the conclusions.

df_sample <- read_csv("../results/samples_hk_study.csv")
caseid_labid <- read_csv("../../2020-09-01_COVID_NGS_pipeline/NGS_data_input/labID_caseID_merged.csv", col_types = cols(.default = "c"))
names(df_sample)[1] <- "lab_ID"

df_trans_pair <- read_csv("../results/identified_transmission_pairs.csv")
df_trans_pair <- df_trans_pair %>% select(pair, case_id) %>% unique()

df_sample <- left_join(df_sample, caseid_labid)
df_master <- read_excel("../../2020-09-01_COVID_NGS_pipeline/pCloud_data/Leo WGS master file.xlsx")

check_trans_pair <- sapply(df_trans_pair$case_id, function(x){
	case_1 <- strsplit(x, " ")[[1]][1]
	case_2 <- strsplit(x, " ")[[1]][2]
	whp_1 <- df_sample$lab_ID[df_sample$case_ID==case_1]
	whp_2 <- df_sample$lab_ID[df_sample$case_ID==case_2]
	df_tmp <- df_master %>% filter(`Lab ID` %in% c(whp_1, whp_2))
	as.character(df_tmp$`Full Genome Date`)
})
check_trans_pair
check_trans_pair <- apply(check_trans_pair, 2, function(x){x[1]==x[2]})
sum(check_trans_pair)
sum(check_trans_pair)/length(check_trans_pair)

df_control <- read_csv("../results/duplicated_cases.csv")
sapply(unique(df_control$case_id), function(x){
	df_tmp <- df_control %>% filter(case_id == x)
	df_t <- df_master %>% filter(`Lab ID` %in% df_tmp$Sample)
	print(x)
	print(df_tmp$same_sample)
	print(as.character(df_t$`NovaSeq Date`))
	return("")
})

