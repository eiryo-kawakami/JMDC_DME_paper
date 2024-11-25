library(lme4)
library(readr)
library(dplyr)
library(investr)
library(uwot)
library(amap)


cols <- c("#601986","#7F4F21","#D70051","#009A3E")

surv_func_test <- read_tsv(paste0("RSF_models_maculopathy_risk_score_mean.csv"))

# surv_func_test$crisis_1year <- 1 - (surv_func_test %>% pull("12.0"))
# surv_func_test$crisis_5year <- 1 - (surv_func_test %>% pull("60.0"))
surv_func_test$month_outcome <- -surv_func_test$month_maculopathy

surv_func_test_DME <- surv_func_test %>% filter(is_mac==1)
surv_func_test_nonDME <- surv_func_test %>% filter(is_mac==0)

logisticNLS_DME <- nls(risk_score ~ SSlogis(month_outcome, A, B, C), data= surv_func_test_DME)

newdata <- data.frame(month_outcome = seq(-177, 0, by = 1))

pred_DME <- as_tibble(predFit(object=logisticNLS_DME, newdata=newdata, interval="confidence")) %>% mutate(month_outcome = newdata$month_outcome)

logisticNLS_nonDME <- nls(risk_score ~ SSlogis(month_outcome, A, B, C), data= surv_func_test_nonDME)

pred_nonDME <- as_tibble(predFit(object=logisticNLS_nonDME, newdata=newdata, interval="confidence")) %>% mutate(month_outcome = newdata$month_outcome)

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

alph = 0.4

## plot test data ##
p <- ggplot(data=surv_func_test%>%filter(is_mac==1))
# p <- p + geom_ribbon(aes(x=month_outcome,ymin=lwr,ymax=upr),fill="#EF857A",data=interval_DME,alpha=0.2)
p <- p + geom_line(aes(x=month_outcome,y=risk_score,group=`加入者id`),colour="#EF857A", alpha=alph, linewidth=1)
p <- p + geom_point(aes(x=month_outcome,y=risk_score),colour="#EF857A", alpha=alph, shape=1, size=1, stroke=1)
p <- p + geom_smooth(aes(month_outcome,y=risk_score),colour="black", fill="#EF857A")
p <- p + theme_tmp
p <- p + scale_y_continuous(limits=c(0,310))
p <- p + scale_x_continuous(limits=c(-168,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-168, 0, 12))

ggsave(paste0("RSF_maculopathy_risk_score_test_DME.pdf"),plot=p,height=3,width=6)

p <- ggplot(data=surv_func_test%>%filter(is_mac==0))
# p <- p + geom_ribbon(aes(x=month_outcome,ymin=lwr,ymax=upr),fill="#82ACC9",data=interval_nonDME,alpha=0.2)
p <- p + geom_line(aes(x=month_outcome,y=risk_score,group=`加入者id`),colour="#82ACC9", alpha=alph, linewidth=1)
p <- p + geom_point(aes(x=month_outcome,y=risk_score),colour="#82ACC9", alpha=alph, shape=1, size=1,stroke=1)
p <- p + geom_smooth(aes(month_outcome,y=risk_score),colour="black", fill="#82ACC9")
p <- p + theme_tmp
p <- p + scale_y_continuous(limits=c(0,310))
p <- p + scale_x_continuous(limits=c(-168,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-168, 0, 12))

ggsave(paste0("RSF_maculopathy_risk_score_test_nonDME.pdf"),plot=p,height=3,width=6)

nlmer_DME <- nlmer(risk_score/300 ~ SSlogis(month_outcome,A, B, C) ~ (A | `加入者id`) + (B | `加入者id`) + (C | `加入者id`), surv_func_test_DME, start = c(A =coef(logisticNLS_DME)["A"][[1]]/300, B = coef(logisticNLS_DME)["B"][[1]], C = coef(logisticNLS_DME)["C"][[1]]))

nlmer_nonDME <- nlmer(risk_score/300 ~ SSlogis(month_outcome,A, B, C) ~ (A | `加入者id`) + (B | `加入者id`) + (C | `加入者id`), surv_func_test_nonDME, start = c(A =coef(logisticNLS_nonDME)["A"][[1]]/300, B = coef(logisticNLS_nonDME)["B"][[1]], C = coef(logisticNLS_nonDME)["C"][[1]]))

SSlogis_parameter <- ranef(nlmer_DME)$`加入者id`
d <- Dist(SSlogis_parameter,method="euclidean")
hc <- hclust(d,method="ward.D2")
dd.col <- as.dendrogram(hc)

cluster_patient <- tibble(Patient_ID=rownames(SSlogis_parameter),cluster=cutree(hc,3))

cluster_patient %>% write_tsv("RSF_maculopathy_risk_score_DME_cluster.txt")

for (i in 1:3){
	patient_selected <- cluster_patient %>% filter(cluster==i) %>% pull(Patient_ID)
	print(length(patient_selected))
	p <- ggplot(data=surv_func_test_DME%>%filter(`加入者id` %in% patient_selected))
	p <- p + geom_line(aes(x=month_outcome,y=risk_score,group=`加入者id`),colour=cols[i], alpha=alph, linewidth=1)
	p <- p + geom_point(aes(x=month_outcome,y=risk_score),colour=cols[i], alpha=alph, shape=1, size=1,stroke=1)
	p <- p + geom_smooth(aes(month_outcome,y=risk_score),colour="black", fill=cols[i])
	p <- p + theme_tmp
	p <- p + scale_y_continuous(limits=c(0,310))
	p <- p + scale_x_continuous(limits=c(-168,0), breaks = seq(-120, 0, 60), minor_breaks = seq(-168, 0, 12))

	ggsave(paste0("RSF_maculopathy_risk_score_test_cluster",i,".pdf"),plot=p,height=3,width=6)
}
