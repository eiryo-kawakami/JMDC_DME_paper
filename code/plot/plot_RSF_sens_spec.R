library(ggplot2)
library(readr)
library(dplyr)
library(ROCit)
library(caret)

n_top_features <- 20

risk_score <- read_tsv("RSF_models_maculopathy_risk_score_mean.csv")

# surv_func_test_DME <- surv_func_test %>% filter(is_mac==1)
# surv_func_test_nonDME <- surv_func_test %>% filter(is_mac==0)

roc_empirical <- rocit(score = risk_score$risk_score, class = risk_score$is_mac) 
auc_ci <- ciAUC(roc_empirical)

AUC_summary <- c()
measure_summary <- c()

threshold_list <- seq(150,250,10)

for (i in seq(0,120,12)){
	risk_score_interval <- risk_score %>% filter(month_maculopathy >= i)

	roc_empirical <- rocit(score = risk_score_interval$risk_score, class = risk_score_interval$is_mac) 
	auc_ci <- ciAUC(roc_empirical)
	
	AUC_summary <- rbind(AUC_summary, c(-i, auc_ci$AUC, auc_ci$upper,auc_ci$lower))

	n_DME <- nrow(risk_score_interval %>% filter(is_mac==1))
	n_nonDME <- nrow(risk_score_interval %>% filter(is_mac==0))

	for (th in threshold_list){

		risk_score_interval_count <- risk_score_interval %>% filter(risk_score>th) %>% group_by(is_mac) %>% summarise(count = n())

		sens <- ifelse(1 %in% risk_score_interval_count$is_mac, risk_score_interval_count %>% filter(is_mac==1) %>% pull(count) / n_DME, 0)
		spec <- ifelse(0 %in% risk_score_interval_count$is_mac, 1- (risk_score_interval_count %>% filter(is_mac==0) %>% pull(count) / n_nonDME), 0)

		measure_summary <- rbind(measure_summary, c(-i, th, sens, spec))
	}
}

measure_summary <- data.frame(measure_summary)
colnames(measure_summary) <- c("time","threshold","sens","spec")
measure_summary$threshold <- factor(measure_summary$threshold)
measure_summary %>% write_tsv("RSF_maculopathy_risk_score_sens_spec.txt")

AUC_summary <- data.frame(AUC_summary)
colnames(AUC_summary) <- c("time","AUC","AUC_upper","AUC_lower")

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
p <- p + geom_ribbon(aes(x=time,ymin=AUC_lower,ymax=AUC_upper),fill="black",alpha=0.3,data=AUC_summary)
p <- p + geom_line(aes(x=time,y=AUC),colour="black",linewidth=1,data=AUC_summary)
p <- p + geom_point(aes(x=time,y=AUC),shape=21,colour="black",fill="white",stroke=1,data=AUC_summary)
p <- p + geom_line(aes(x=time,y=sens,linetype=threshold),colour="#EF857A",linewidth=1,data=measure_summary %>% filter(threshold %in% c("150","170","190")))
p <- p + geom_point(aes(x=time,y=sens),shape=21,colour="#EF857A",fill="white",stroke=1,data=measure_summary %>% filter(threshold %in% c("150","170","190")))
p <- p + geom_line(aes(x=time,y=spec,linetype=threshold),colour="#82ACC9",linewidth=1,data=measure_summary %>% filter(threshold %in% c("150","170","190")))
p <- p + geom_point(aes(x=time,y=spec),shape=21,colour="#82ACC9",fill="white",stroke=1,data=measure_summary %>% filter(threshold %in% c("150","170","190")))
p <- p + scale_x_continuous(limits=c(NA,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-126, 0, 6))
p <- p + scale_y_continuous(limits=c(0,1), breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.1))
p <- p + theme_tmp

ggsave(paste0("RSF_maculopathy_risk_score_performance.pdf"),plot=p,height=2.5,width=6)

