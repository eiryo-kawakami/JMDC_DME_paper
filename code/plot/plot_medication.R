library(ggplot2)
library(readr)
library(dplyr)
library(paletteer)


cols <- c("#4E79A7FF","#F28E2BFF","#E15759FF","#76B7B2FF","#59A14FFF","#EDC948FF","#B07AA1FF","#FF9DA7FF","#9C755FFF")

surv_func_test <- read_tsv("RSF_maculopathy_risk_survival_func_mean_test.txt")
nonDME_test <- surv_func_test %>% filter(is_mac==0) %>% pull(`加入者id`) %>% unique()

DME_cluster <- read_tsv("RSF_maculopathy_risk_score_DME_cluster.txt")

data_all <- read_csv("before_interpolate.csv")

data_DME_test <- data_all %>% filter(`加入者id` %in% DME_cluster$Patient_ID)
data_DME_test_cluster <- left_join(data_DME_test,DME_cluster,by=c("加入者id"="Patient_ID"))

data_nonDME_test_cluster <- data_all %>% filter(`加入者id` %in% nonDME_test)
data_nonDME_test_cluster$cluster <- "nonDME"

data_test_cluster <- rbind(data_DME_test_cluster,data_nonDME_test_cluster)

data_test_cluster$cluster <- factor(data_test_cluster$cluster,levels=c("3","1","2","nonDME"))
data_test_cluster$month_outcome <- -data_test_cluster$month_maculopathy

medication_data <- read_tsv("iyakuhin_mac_ctl_month.txt")
medication_category <- read_csv("糖尿病薬リスト_医薬品コード_SGLT2追加.csv")
medication_group <- unique(medication_category$group)

medication_data <- left_join(medication_data,DME_cluster,by=c("加入者id"="Patient_ID"))
medication_data$cluster[is.na(medication_data$cluster)] <- "nonDME"
medication_data$cluster <- factor(medication_data$cluster,levels=c("3","1","2","nonDME"))

medication_test <- c()

for (i in seq(12,120,12)){
	medication_data_interval <- medication_data %>% filter((month_maculopathy >= (i-12)) & (month_maculopathy < i))
	data_test_cluster_interval <- data_test_cluster %>% filter((month_maculopathy >= (i-12)) & (month_maculopathy < i))

	for (m in medication_group){
		medication_list <- medication_category %>% filter(group==m) %>% pull("medicine")
		for (c in c("3","2","1","nonDME")){
			all_patient <- data_test_cluster_interval %>% filter(cluster==c) %>% pull("加入者id") %>% unique()
			freq <- medication_data_interval %>% filter(cluster==c) %>% filter(医薬品名 %in% medication_list) %>% filter(加入者id %in% all_patient) %>% pull("加入者id") %>% unique() %>% length() / length(all_patient)

			medication_test <- rbind(medication_test,c(freq,m,-i+6,c))
		}
	}
}

medication_data_interval <- medication_data %>% filter(month_maculopathy >= 120)
data_test_cluster_interval <- data_test_cluster %>% filter(month_maculopathy >= 120)

for (m in medication_group){
	medication_list <- medication_category %>% filter(group==m) %>% pull("medicine")
	for (c in c("3","1","2","nonDME")){
		all_patient <- data_test_cluster_interval %>% filter(cluster==c) %>% pull("加入者id") %>% unique()
		freq <- medication_data_interval %>% filter(cluster==c) %>% filter(医薬品名 %in% medication_list) %>% filter(加入者id %in% all_patient) %>% pull("加入者id") %>% unique() %>% length() / length(all_patient)

		medication_test <- rbind(medication_test,c(freq,m,-126,c))
	}
}

medication_test <- data.frame(medication_test)
colnames(medication_test) <- c("freq","medication","time","cluster")

medication_test$cluster <- factor(medication_test$cluster,levels=c("3","1","2","nonDME"))
medication_test$time <- as.numeric(medication_test$time)
medication_test$freq <- as.numeric(medication_test$freq)
medication_test$medication <- factor(medication_test$medication,levels=unique(medication_test$medication))


theme_tmp <- theme(
  axis.title = element_text(family = "sans", size = 16), 
  axis.text = element_text(family = "sans", size = 16), 
  legend.title = element_text(family = "sans", size = 16), 
  legend.text = element_text(family = "sans", size = 16), 
  panel.background = element_blank(), 
  panel.border = element_rect(fill = NA, color = "black"), 
  panel.grid.major = element_line(color = "black", linewidth = 0.2), 
  panel.grid.minor = element_line(color = "black", linetype = 2, linewidth = 0.2)
)

p <- ggplot(data=medication_test)
p <- p + geom_line(aes(x=time,y=freq,colour=medication),linewidth=1)
p <- p + geom_point(aes(x=time,y=freq,colour=medication),shape=21,fill="white",stroke=1)
p <- p + scale_colour_manual(values=cols)
p <- p + scale_x_continuous(limits=c(NA,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-120, 0, 6))
p <- p + scale_y_continuous(limits=c(0,1), breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.1))
p <- p + theme_tmp
p <- p + facet_wrap(~ cluster, ncol=2)

ggsave(paste0("RSF_maculopathy_medication_risk_score_cluster.pdf"),plot=p,height=3.5,width=8)
