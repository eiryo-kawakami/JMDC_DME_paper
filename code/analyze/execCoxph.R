require(survival)
require(coxme)
library(readr)
library(dplyr)


data_all <- read_csv("before_interpolate.csv")

binary_vars <- c("既往歴","自覚症状","他覚症状","心電図(所見の有無)","喫煙","食べ方2(就寝前)","食べ方3(夜食/間食)","食習慣","睡眠","服薬1(血圧)","服薬2(血糖)","服薬3(脂質)","既往歴1(脳血管)","既往歴2(心血管)","既往歴3(腎不全･人工透析)","貧血","20歳からの体重変化","30分以上の運動習慣","歩行又は身体活動","歩行速度","1年間の体重変化","保健指導の希望")

res_summary <- c()

for (i in 2:(ncol(data_all)-2)){

	print(colnames(data_all)[i])

	data_ind <- data.frame(data_all[,c("加入者id","is_mac","month_maculopathy")],value=data_all[,i])
	colnames(data_ind) <- c("patient_ID","is_mac","month_maculopathy","value")

	if (colnames(data_all)[i] %in% binary_vars){
		data_ind$value <- 2 - data_ind$value
	}

	if (colnames(data_all)[i] == "飲酒"){
		data_ind$value <- 3 - data_ind$value
	}

	nonzero_D <- sum(data_ind[data_ind$is_mac==1,"value"]!=0, na.rm=T)
	nonzero_ND <- sum(data_ind[data_ind$is_mac==0,"value"]!=0, na.rm=T)

	res <- coxme(Surv(month_maculopathy, is_mac) ~ value + (1|patient_ID), data = data_ind)
	res_summary <- rbind(res_summary,c(colnames(data_all)[i],summary(res)$coefficients["value",c("coef","z","p")]))
	# } else {
	# 	res_summary <- rbind(res_summary,c(colnames(data_all)[i],NA,NA,NA))
	# }
}

colnames(res_summary) <- c("variable","coef","z_score","p_value")
res_summary <- data.frame(res_summary)
res_summary$FDR <- p.adjust(res_summary$p_value)

write.table(res_summary,file="JMDC_DME_Coxme.txt",sep="\t",quote=FALSE,row.names=FALSE)


res_summary <- c()

for (i in 3:(ncol(data_all)-2)){

	print(colnames(data_all)[i])

	data_ind <- data.frame(data_all[,c("加入者id","is_mac","month_maculopathy","age")],value=data_all[,i])
	colnames(data_ind) <- c("patient_ID","is_mac","month_maculopathy","age","value")

	if (colnames(data_all)[i] %in% binary_vars){
		data_ind$value <- 2 - data_ind$value
	}

	if (colnames(data_all)[i] == "飲酒"){
		data_ind$value <- 3 - data_ind$value
	}

	nonzero_D <- sum(data_ind[data_ind$is_mac==1,"value"]!=0, na.rm=T)
	nonzero_ND <- sum(data_ind[data_ind$is_mac==0,"value"]!=0, na.rm=T)

	res <- coxme(Surv(month_maculopathy, is_mac) ~ value + age + (1|patient_ID), data = data_ind)
	res_summary <- rbind(res_summary,c(colnames(data_all)[i],summary(res)$coefficients["value",c("coef","z","p")]))
	# } else {
	# 	res_summary <- rbind(res_summary,c(colnames(data_all)[i],NA,NA,NA))
	# }
}

colnames(res_summary) <- c("variable","coef","z_score","p_value")
res_summary <- data.frame(res_summary)
res_summary$FDR <- p.adjust(res_summary$p_value)

write.table(res_summary,file="JMDC_DME_Coxme_age_adjusted.txt",sep="\t",quote=FALSE,row.names=FALSE)
