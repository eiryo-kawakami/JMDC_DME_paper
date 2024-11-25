library(readr)
library(ggplot2)
library(dplyr)


varimp_all <- read_tsv("RSF_maculopathy_risk_varimp_descending.txt")
varimp_cluster1 <- read_tsv("RSF_maculopathy_risk_varimp_cluster1_descending.txt")
varimp_cluster2 <- read_tsv("RSF_maculopathy_risk_varimp_cluster2_descending.txt")
varimp_cluster3 <- read_tsv("RSF_maculopathy_risk_varimp_cluster3_descending.txt")

varimp_all_top10 <- varimp_all %>% top_n(n=10, median)
varimp_all_top10$index <- factor(varimp_all_top10$index,levels=rev(varimp_all_top10$index))

varimp_cluster1_top10 <- varimp_cluster1 %>% top_n(n=10, median)
varimp_cluster1_top10$index <- factor(varimp_cluster1_top10$index,levels=rev(varimp_cluster1_top10$index))

varimp_cluster2_top10 <- varimp_cluster2 %>% top_n(n=10, median)
varimp_cluster2_top10$index <- factor(varimp_cluster2_top10$index,levels=rev(varimp_cluster2_top10$index))

varimp_cluster3_top10 <- varimp_cluster3 %>% top_n(n=10, median)
varimp_cluster3_top10$index <- factor(varimp_cluster3_top10$index,levels=rev(varimp_cluster3_top10$index))

syobyo_master <- read_csv("傷病マスタ.csv")
disease <- syobyo_master$`標準病名`

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

varimp_all_top10$category <- factor(ifelse(varimp_all_top10$index %in% disease, "disease", "checkup"),levels=c("checkup","disease"))
varimp_cluster1_top10$category <- ifelse(varimp_cluster1_top10$index %in% disease, "disease", "checkup")
varimp_cluster2_top10$category <- ifelse(varimp_cluster2_top10$index %in% disease, "disease", "checkup")
varimp_cluster3_top10$category <- ifelse(varimp_cluster3_top10$index %in% disease, "disease", "checkup")

p <- ggplot(data=varimp_all_top10)
p <- p + geom_bar(aes(x=index,y=median,fill=category),width=0.8,position = position_dodge(width = 0.8),stat="identity",colour="black")
p <- p + scale_fill_manual(values=c("white","black"))
p <- p + theme_tmp
p <- p + coord_flip()
p <- p + theme(aspect.ratio = 3)

ggsave("RSF_maculopathy_risk_varimp.pdf",plot=p,height=5,width=6)

p <- ggplot(data=varimp_cluster1_top10)
p <- p + geom_bar(aes(x=index,y=median,fill=category),width=0.8,position = position_dodge(width = 0.8),stat="identity",colour="black")
p <- p + scale_fill_manual(values=c("white","black"))
p <- p + theme_tmp
p <- p + coord_flip()
p <- p + theme(aspect.ratio = 3)

ggsave("RSF_maculopathy_risk_varimp_cluster1.pdf",plot=p,height=5,width=6)

p <- ggplot(data=varimp_cluster2_top10)
p <- p + geom_bar(aes(x=index,y=median,fill=category),width=0.8,position = position_dodge(width = 0.8),stat="identity",colour="black")
p <- p + scale_fill_manual(values=c("white","black"))
p <- p + theme_tmp
p <- p + coord_flip()
p <- p + theme(aspect.ratio = 3)

ggsave("RSF_maculopathy_risk_varimp_cluster2.pdf",plot=p,height=5,width=6)

p <- ggplot(data=varimp_cluster3_top10)
p <- p + geom_bar(aes(x=index,y=median,fill=category),width=0.8,position = position_dodge(width = 0.8),stat="identity",colour="black")
p <- p + scale_fill_manual(values=c("white","black"))
p <- p + theme_tmp
p <- p + coord_flip()
p <- p + theme(aspect.ratio = 3)

ggsave("RSF_maculopathy_risk_varimp_cluster3.pdf",plot=p,height=5,width=6)
