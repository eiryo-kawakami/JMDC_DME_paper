library(ggplot2)
library(readr)
library(dplyr)
library(rstatix)


cols <- c("#601986","#D70051","#7F4F21","#82ACC9")

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

p <- ggplot()
p <- p + geom_smooth(data=data_test_cluster, aes(x=month_outcome, y=hba1c, colour=cluster, fill=cluster))
p <- p + scale_colour_manual(values=cols[c(2,1,3,4)])
p <- p + scale_fill_manual(values=cols[c(2,1,3,4)])
p <- p + theme_tmp
p <- p + scale_x_continuous(limits=c(-168,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-168, 0, 12))

ggsave(paste0("RSF_maculopathy_HbA1c_risk_score_cluster_test_DME.pdf"),plot=p,height=3,width=6)

p <- ggplot()
p <- p + geom_smooth(data=data_test_cluster, aes(x=month_outcome, y=空腹時血糖, colour=cluster, fill=cluster))
p <- p + scale_colour_manual(values=cols[c(2,1,3,4)])
p <- p + scale_fill_manual(values=cols[c(2,1,3,4)])
p <- p + theme_tmp
p <- p + scale_x_continuous(limits=c(-168,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-168, 0, 12))

ggsave(paste0("RSF_maculopathy_FBG_risk_score_cluster_test_DME.pdf"),plot=p,height=3,width=6)

p <- ggplot()
p <- p + geom_smooth(data=data_test_cluster, aes(x=month_outcome, y=`γ-gt(γ-gtp)`, colour=cluster, fill=cluster))
p <- p + scale_colour_manual(values=cols[c(2,1,3,4)])
p <- p + scale_fill_manual(values=cols[c(2,1,3,4)])
p <- p + theme_tmp
p <- p + scale_x_continuous(limits=c(-168,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-168, 0, 12))

ggsave(paste0("RSF_maculopathy_GGT_risk_score_cluster_test_DME.pdf"),plot=p,height=3,width=6)

p <- ggplot()
p <- p + geom_smooth(data=data_test_cluster, aes(x=month_outcome, y=`収縮期血圧`, colour=cluster, fill=cluster))
p <- p + scale_colour_manual(values=cols[c(2,1,3,4)])
p <- p + scale_fill_manual(values=cols[c(2,1,3,4)])
p <- p + theme_tmp
p <- p + scale_x_continuous(limits=c(-168,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-168, 0, 12))

ggsave(paste0("RSF_maculopathy_SBP_risk_score_cluster_test_DME.pdf"),plot=p,height=3,width=6)

p <- ggplot()
p <- p + geom_smooth(data=data_test_cluster, aes(x=month_outcome, y=bmi, colour=cluster, fill=cluster))
p <- p + scale_colour_manual(values=cols[c(2,1,3,4)])
p <- p + scale_fill_manual(values=cols[c(2,1,3,4)])
p <- p + theme_tmp
p <- p + scale_x_continuous(limits=c(-168,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-168, 0, 12))

ggsave(paste0("RSF_maculopathy_BMI_risk_score_cluster_test_DME.pdf"),plot=p,height=3,width=6)

p <- ggplot()
p <- p + geom_smooth(data=data_test_cluster, aes(x=month_outcome, y=`血色素量(ヘモグロビン値)`, colour=cluster, fill=cluster))
p <- p + scale_colour_manual(values=cols[c(2,1,3,4)])
p <- p + scale_fill_manual(values=cols[c(2,1,3,4)])
p <- p + theme_tmp
p <- p + scale_x_continuous(limits=c(-168,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-168, 0, 12))

ggsave(paste0("RSF_maculopathy_Hemoglobin_risk_score_cluster_test_DME.pdf"),plot=p,height=3,width=6)

p <- ggplot()
p <- p + geom_smooth(data=data_test_cluster, aes(x=month_outcome, y=`ヘマトクリット値`, colour=cluster, fill=cluster))
p <- p + scale_colour_manual(values=cols[c(2,1,3,4)])
p <- p + scale_fill_manual(values=cols[c(2,1,3,4)])
p <- p + theme_tmp
p <- p + scale_x_continuous(limits=c(-168,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-168, 0, 12))

ggsave(paste0("RSF_maculopathy_Hematocrit_risk_score_cluster_test_DME.pdf"),plot=p,height=3,width=6)

nonDME_cluster <- data.frame(Patient_ID=nonDME_test,cluster="nonDME")
test_cluster <- rbind(DME_cluster,nonDME_cluster)

baseline_data <- read_csv("matched_list.csv")
baseline_cluster <- left_join(test_cluster,baseline_data,by=c("Patient_ID"="pid"))
baseline_cluster$cluster <- factor(baseline_cluster$cluster,levels=c("3","1","2","nonDME"))

stat.test <- baseline_cluster %>% t_test(age ~ cluster) %>% filter(p.adj.signif!="ns")

p <- ggplot()
p <- p + geom_boxplot(data=baseline_cluster, aes(x=cluster, y=age, colour=cluster))
p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.03, y.position = 80, step.increase = 0.1)
p <- p + scale_colour_manual(values=cols[c(2,1,3,4)])
# p <- p + scale_fill_manual(values=cols[c(2,1,3,4)])
p <- p + theme_bw(base_size=16)

ggsave(paste0("RSF_maculopathy_age_at_start_risk_score_cluster_test_DME.pdf"),plot=p,height=3.5,width=4)

baseline_cluster$age_diabetes <- baseline_cluster$age - round((baseline_cluster$durling_DM - baseline_cluster$period)/12)

stat.test <- baseline_cluster %>% t_test(age_diabetes ~ cluster) %>% filter(p.adj.signif!="ns")

p <- ggplot()
p <- p + geom_boxplot(data=baseline_cluster, aes(x=cluster, y=age, colour=cluster))
p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.03, y.position = 80, step.increase = 0.1)
p <- p + scale_colour_manual(values=cols[c(2,1,3,4)])
# p <- p + scale_fill_manual(values=cols[c(2,1,3,4)])
p <- p + theme_bw(base_size=16)

ggsave(paste0("RSF_maculopathy_age_at_start_risk_score_cluster_test_DME.pdf"),plot=p,height=3.5,width=4)

stat.test <- baseline_cluster %>% t_test(durling_DM ~ cluster) %>% filter(p.adj.signif!="ns")

p <- ggplot()
p <- p + geom_boxplot(data=baseline_cluster, aes(x=cluster, y=durling_DM, colour=cluster))
p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.03, y.position = 400, step.increase = 0.1)
p <- p + scale_colour_manual(values=cols[c(2,1,3,4)])
# p <- p + scale_fill_manual(values=cols[c(2,1,3,4)])
p <- p + theme_bw(base_size=16)

ggsave(paste0("RSF_maculopathy_durationDM_risk_score_cluster_test_DME.pdf"),plot=p,height=3.5,width=4)

cluster_rate_M <- baseline_cluster %>% filter(sex=="M") %>% filter(cluster!="nonDME") %>% group_by(cluster) %>% summarise(N=n())
cluster_rate_M$cluster <- factor(cluster_rate_M$cluster,levels=c("3","1","2"))
cluster_rate_M <- cluster_rate_M[order(cluster_rate_M$cluster),]

p <- ggplot(cluster_rate_M)
p <- p + geom_bar(aes(x="", y=N, fill=cluster),stat="identity", width=1, color="white")
p <- p + coord_polar("y", start=0)
p <- p + theme_void()
p <- p + theme(legend.position="none")
p <- p + scale_fill_manual(values=cols[c(3,1,2)])

ggsave(paste0("RSF_maculopathy_genderM_risk_score_cluster_test_DME.pdf"),plot=p,height=3,width=4)


cluster_rate_F <- baseline_cluster %>% filter(sex=="F") %>% filter(cluster!="nonDME") %>% group_by(cluster) %>% summarise(N=n())
cluster_rate_F$cluster <- factor(cluster_rate_F$cluster,levels=c("3","1","2"))
cluster_rate_F <- cluster_rate_F[order(cluster_rate_F$cluster),]

p <- ggplot(cluster_rate_F)
p <- p + geom_bar(aes(x="", y=N, fill=cluster),stat="identity", width=1, color="white")
p <- p + coord_polar("y", start=0)
p <- p + theme_void()
p <- p + theme(legend.position="none")
p <- p + scale_fill_manual(values=cols[c(3,1,2)])

ggsave(paste0("RSF_maculopathy_genderF_risk_score_cluster_test_DME.pdf"),plot=p,height=3,width=4)


disease_list <- c("糖尿病性腎症","遠視性乱視","糖尿病性末梢神経障害")
disease_freq <- c()

for (i in seq(12,120,12)){
	data_test_cluster_interval <- data_test_cluster %>% filter((month_maculopathy >= (i-12)) & (month_maculopathy < i))

	for (d in disease_list){
		disease_freq <- rbind(disease_freq,c(d,-i+6,"2",data_test_cluster_interval %>% filter(cluster=="2") %>% pull(d) %>% sum()/nrow(data_test_cluster_interval %>% filter(cluster=="2"))))
		disease_freq <- rbind(disease_freq,c(d,-i+6,"1",data_test_cluster_interval %>% filter(cluster=="1") %>% pull(d) %>% sum()/nrow(data_test_cluster_interval %>% filter(cluster=="1"))))
		disease_freq <- rbind(disease_freq,c(d,-i+6,"3",data_test_cluster_interval %>% filter(cluster=="3") %>% pull(d) %>% sum()/nrow(data_test_cluster_interval %>% filter(cluster=="3"))))
		disease_freq <- rbind(disease_freq,c(d,-i+6,"nonDME",data_test_cluster_interval %>% filter(cluster=="nonDME") %>% pull(d) %>% sum()/nrow(data_test_cluster_interval %>% filter(cluster=="nonDME"))))
	}
}

data_test_cluster_interval <- data_test_cluster %>% filter(month_maculopathy >= 120)

for (d in disease_list){
	disease_freq <- rbind(disease_freq,c(d,-126,"2",data_test_cluster_interval %>% filter(cluster=="2") %>% pull(d) %>% sum()/nrow(data_test_cluster_interval %>% filter(cluster=="2"))))
	disease_freq <- rbind(disease_freq,c(d,-126,"1",data_test_cluster_interval %>% filter(cluster=="1") %>% pull(d) %>% sum()/nrow(data_test_cluster_interval %>% filter(cluster=="1"))))
	disease_freq <- rbind(disease_freq,c(d,-126,"3",data_test_cluster_interval %>% filter(cluster=="3") %>% pull(d) %>% sum()/nrow(data_test_cluster_interval %>% filter(cluster=="3"))))
	disease_freq <- rbind(disease_freq,c(d,-126,"nonDME",data_test_cluster_interval %>% filter(cluster=="nonDME") %>% pull(d) %>% sum()/nrow(data_test_cluster_interval %>% filter(cluster=="nonDME"))))
}

disease_freq <- data.frame(disease_freq)
colnames(disease_freq) <- c("disease","time","cluster","freq")
disease_freq$cluster <- factor(disease_freq$cluster,levels=c("3","1","2","nonDME"))
disease_freq$time <- as.numeric(disease_freq$time)
disease_freq$freq <- as.numeric(disease_freq$freq)

for (d in disease_list){
	disease_freq_ind <- disease_freq %>% filter(disease==d)

	p <- ggplot(data=disease_freq_ind)
	p <- p + geom_line(aes(x=time,y=freq,colour=cluster),linewidth=1)
	p <- p + geom_point(aes(x=time,y=freq,colour=cluster),shape=21,fill="white",stroke=1)
	p <- p + scale_colour_manual(values=cols[c(2,1,3,4)])
	p <- p + scale_x_continuous(limits=c(NA,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-120, 0, 6))
	# p <- p + scale_y_continuous(limits=c(0,1), breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.1))
	p <- p + theme_tmp

	ggsave(paste0("RSF_maculopathy_",d,"_risk_score_cluster_test_DME.pdf"),plot=p,height=3,width=6)
}


## 尿検査のカウント
urinary_list <- c("尿蛋白(定性)","尿糖")
urinary_test <- c()

for (i in seq(12,120,12)){
	data_test_cluster_interval <- data_test_cluster %>% filter((month_maculopathy >= (i-12)) & (month_maculopathy < i))

	for (u in urinary_list){
		for (c in c("2","1","3","nonDME")){
			tmp_data <- data_test_cluster_interval %>% filter(cluster==c) %>% pull(u) %>% table() %>% as.data.frame()
			colnames(tmp_data) <- c("level","freq")
			tmp_data$freq <- tmp_data$freq / nrow(data_test_cluster_interval %>% filter(!is.na(data_test_cluster_interval[u] %>% unlist())) %>% filter(cluster==c))
			tmp_data$test <- u
			tmp_data$time <- -i+6
			tmp_data$cluster <- c
			tmp_data$level <- as.numeric(tmp_data$level)

			for (l in c(1:5)){
				if (!(l %in% tmp_data$level)){
					tmp_data <- rbind(tmp_data,c(l,0,u,-i+6,c))
				}
			}
			urinary_test <- rbind(urinary_test,tmp_data)
		}
	}
}

data_test_cluster_interval <- data_test_cluster %>% filter(month_maculopathy >= 120)

for (u in urinary_list){
	for (c in c("2","1","3","nonDME")){
		tmp_data <- data_test_cluster_interval %>% filter(cluster==c) %>% pull(u) %>% table() %>% as.data.frame()
		colnames(tmp_data) <- c("level","freq")
		tmp_data$freq <- tmp_data$freq / nrow(data_test_cluster_interval %>% filter(!is.na(data_test_cluster_interval[u] %>% unlist())) %>% filter(cluster==c))
		tmp_data$test <- u
		tmp_data$time <- -126
		tmp_data$cluster <- c
		tmp_data$level <- as.numeric(tmp_data$level)

		for (l in c(1:5)){
			if (!(l %in% tmp_data$level)){
				tmp_data <- rbind(tmp_data,c(l,0,u,-126,c))
			}
		}
		urinary_test <- rbind(urinary_test,tmp_data)
	}
}


urinary_test$cluster <- factor(urinary_test$cluster,levels=c("3","1","2","nonDME"))
urinary_test$time <- as.numeric(urinary_test$time)
urinary_test$freq <- as.numeric(urinary_test$freq)
urinary_test$level <- factor(urinary_test$level)

for (u in urinary_list){
	urinary_test_ind <- urinary_test %>% filter(test==u)

	p <- ggplot(data=urinary_test_ind)
	p <- p + geom_area(aes(x=time,y=freq,fill=level),colour="black", alpha=0.7)
	p <- p + scale_fill_brewer(palette="Oranges")
	p <- p + scale_x_continuous(limits=c(NA,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-120, 0, 6))
	# p <- p + scale_y_continuous(limits=c(0,1), breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.1))
	p <- p + theme_tmp
	p <- p + facet_wrap(~ cluster, ncol=2)

	ggsave(paste0("RSF_maculopathy_",u,"_risk_score_cluster_test_DME.pdf"),plot=p,height=4,width=7)
}


## アルコール摂取のカウント
alcohol <- c()

for (i in seq(12,120,12)){
	data_test_cluster_interval <- data_test_cluster %>% filter((month_maculopathy >= (i-12)) & (month_maculopathy < i))

	for (c in c("2","1","3","nonDME")){
		tmp_data <- data_test_cluster_interval %>% filter(cluster==c) %>% pull("飲酒") %>% table() %>% as.data.frame()
		colnames(tmp_data) <- c("level","freq")
		tmp_data$freq <- tmp_data$freq / nrow(data_test_cluster_interval %>% filter(!is.na(data_test_cluster_interval["飲酒"] %>% unlist())) %>% filter(cluster==c))
		tmp_data$test <- "飲酒"
		tmp_data$time <- -i+6
		tmp_data$cluster <- c
		tmp_data$level <- as.numeric(tmp_data$level)

		for (l in c(1:3)){
			if (!(l %in% tmp_data$level)){
				tmp_data <- rbind(tmp_data,c(l,0,"飲酒",-i+6,c))
			}
		}
		alcohol <- rbind(alcohol,tmp_data)
	}
}

data_test_cluster_interval <- data_test_cluster %>% filter(month_maculopathy >= 120)

for (c in c("2","1","3","nonDME")){
	tmp_data <- data_test_cluster_interval %>% filter(cluster==c) %>% pull("飲酒") %>% table() %>% as.data.frame()
	colnames(tmp_data) <- c("level","freq")
	tmp_data$freq <- tmp_data$freq / nrow(data_test_cluster_interval %>% filter(!is.na(data_test_cluster_interval["飲酒"] %>% unlist())) %>% filter(cluster==c))
	tmp_data$test <- "飲酒"
	tmp_data$time <- -126
	tmp_data$cluster <- c
	tmp_data$level <- as.numeric(tmp_data$level)

	for (l in c(1:3)){
		if (!(l %in% tmp_data$level)){
			tmp_data <- rbind(tmp_data,c(l,0,"飲酒",-126,c))
		}
	}
	alcohol <- rbind(alcohol,tmp_data)
}


alcohol$cluster <- factor(alcohol$cluster,levels=c("3","1","2","nonDME"))
alcohol$time <- as.numeric(alcohol$time)
alcohol$freq <- as.numeric(alcohol$freq)
alcohol$level <- factor(alcohol$level)

p <- ggplot(data=alcohol)
p <- p + geom_area(aes(x=time,y=freq,fill=level),colour="black", alpha=0.7)
p <- p + scale_fill_brewer(palette="Purples")
p <- p + scale_x_continuous(limits=c(NA,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-120, 0, 6))
# p <- p + scale_y_continuous(limits=c(0,1), breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.1))
p <- p + theme_tmp
p <- p + facet_wrap(~ cluster, ncol=2)

ggsave(paste0("RSF_maculopathy_alcohol_freq_risk_score_cluster_test_DME.pdf"),plot=p,height=4,width=7)
