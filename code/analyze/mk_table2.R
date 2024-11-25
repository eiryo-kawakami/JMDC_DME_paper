library(table1)
library(readr)
library(dplyr)

df <- read_csv("before_interpolate.csv")
df$is_mac <- factor(df$is_mac)

health_checkup_vars <- colnames(df)[c(2:46)]

use_vars <- c()

for (v in health_checkup_vars){
	unique_val <- length(unique(df %>% pull(v)))
	print(paste0(v,(": "),unique_val))
	if (unique_val < 10){
		df[v] <- factor(df %>% pull(v))
	}
}

pvalue <- function(x, ...) {
    # Construct vectors of data y, and groups (strata) g
    y <- unlist(x)
    g <- factor(rep(1:length(x), times=sapply(x, length)))
    if (is.numeric(y)) {
        # For numeric variables, perform a standard 2-sample t-test
        p <- t.test(y ~ g)$p.value
    } else {
        # For categorical variables, perform a chi-squared test of independence
        p <- chisq.test(table(y, g))$p.value
    }
    # Format the p-value, using an HTML entity for the less-than sign.
    # The initial empty string places the output on the line below the variable label.
    c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

table1(~ age + gender + bmi + `腹囲` + `既往歴` + `自覚症状` + `他覚症状` + `収縮期血圧` + `拡張期血圧` + `中性脂肪(トリグリセリド)` + `hdlコレステロール` + `ldlコレステロール` + `got(ast)` + `gpt(alt)` + `γ-gt(γ-gtp)` + `空腹時血糖` + `hba1c` + `尿糖` + `尿蛋白(定性)` + `ヘマトクリット値` + `血色素量(ヘモグロビン値)` + `赤血球数` + `心電図(所見の有無)` + `喫煙` + `食べ方1(早食い等)` + `食べ方2(就寝前)` + `食べ方3(夜食/間食)` + `食習慣` + `飲酒` + `飲酒量` + `睡眠` + `服薬1(血圧)` + `服薬2(血糖)` + `服薬3(脂質)` + `既往歴1(脳血管)` + `既往歴2(心血管)` + `既往歴3(腎不全･人工透析)` + `貧血` + `20歳からの体重変化` + `30分以上の運動習慣` + `歩行又は身体活動` + `歩行速度` + `1年間の体重変化` + `生活習慣の改善` + `保健指導の希望` | is_mac, data=df, overall=F)

