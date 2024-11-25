library(readr)
library(dplyr)


syobyo_master <- read_csv("傷病マスタ.csv")
syobyo_master$variable <- syobyo_master$`標準病名`

data_all <- read_csv("before_interpolate.csv")

disease_vars <- colnames(data_all)[c(47:(ncol(data_all)-2))]

data_disease <- data_all %>% select(all_of(c("is_mac",disease_vars)))
data_disease$is_mac <- factor(data_disease$is_mac)

disease_incidence <- data_disease %>% group_by(is_mac) %>% summarise_all(list(sum)) %>% t() %>% as.data.frame()
disease_incidence <- disease_incidence[c(2:nrow(disease_incidence)),]
colnames(disease_incidence) <- c("Non-DME","DME")
disease_incidence$variable <- rownames(disease_incidence)

disease_incidence <- left_join(disease_incidence,syobyo_master %>% select(all_of(c("variable","icd10大分類コード","icd10小分類コード"))),by="variable")

disease_incidence <- disease_incidence %>% select(c("variable","icd10大分類コード","icd10小分類コード","Non-DME","DME"))

cox_res <- read_tsv("JMDC_DME_Coxme_age_adjusted.txt")

disease_summary <- left_join(disease_incidence,cox_res,by="variable")
disease_summary %>% write_tsv("JMDC_DME_Coxme_age_adjusted_disease_summary.txt")

disease_summary %>% filter(FDR < 0.05) %>% group_by(`icd10大分類コード`) %>% summarise(n = n())
