my_packages <- c("devtools", "readxl", "EpiFunc", "janitor", "lubridate", "skimr", 
                 "gtsummary", "ggplot2", "naniar", "openxlsx", 
                 "multcomp", "tidyverse", "survival", "mfp") 
lapply(my_packages, require, character.only = TRUE)

# load the data 
rm(list=ls())
prs <- read.csv("result/prs.csv", header  = TRUE)

# put data into 7 groups according to radarID
information <- read.csv("data/subgroup/NS_clin_lab_med_phenotypes.csv")

names <- prs$name
information <- information %>% 
  filter(patient_id %in% names)

colnames(prs)[colnames(prs) == "name"] <- "patient_id"

dta <- merge(prs, information, by = "patient_id")

dta <- dta[, -2]

dta <- dta %>%
  mutate(rough_NSgroup = recode(rough_NSgroup, 
                                "SRNS genetic" = "Monogenic SRNS",
                                "SRNS non-genetic" = "Non Monogenic SRNS"))

# compare between subgroups

dta <- dta[dta$rough_NSgroup != "Unclear gen/non - exclude in subgroup analyses",]
dta <- dta[!is.na(dta$first_diagnosis_age), ]
dta$rough_NSgroup <- as.factor(dta$rough_NSgroup)

grps <- c("SSNS", "Monogenic SRNS", "Non Monogenic SRNS", "CFD")
four_nsgrps <- dta[dta$rough_NSgroup %in% grps,]

results <- aggregate(dta$prs_multi_ancestry, list(dta$rough_NSgroup), FUN=summary)
write.csv(results, "result/PRS by NS subgroup summary.csv", row.names = FALSE)

a_res <- aov(prs_multi_ancestry ~ rough_NSgroup, data = dta) 
a_NS_res <- summary(a_res)

# Which exact groups differ?
# Tukey HSD test:
post_test <- glht(a_res,linfct = mcp(rough_NSgroup = "Tukey"))
post_NS_res <- summary(post_test)

summary(post_test)

shpal <- c("navyblue", "orange", "green", "red", "purple", "brown", "cyan")
shpal2 <- c("navyblue", "green", "red", "purple")



prs_boxplot <- function(df, y_var, title_suffix, color_palette) {
  p <- ggplot(df, aes(x = rough_NSgroup, y = !!sym(y_var), fill = rough_NSgroup)) + 
    geom_boxplot(width = 0.5) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10), 
      axis.text.y = element_text(face = "bold"), 
      axis.text.x = element_text(face = "bold", angle = 0), 
      axis.title.x = element_blank(), 
      axis.title.y = element_text(face = "bold"), 
      legend.title = element_blank(),
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major.x = element_blank(), # 隐藏纵向的网格线
      panel.grid.minor.x = element_blank(), # 隐藏次纵向的网格线
      panel.grid.major.y = element_line(color = "grey80"), # 添加水平网格线
      panel.grid.minor.y = element_blank(), # 隐藏次水平网格线
      axis.line = element_line(color = "black") # 显示轴线
    ) +
    labs(y = "Polygenic Risk Score", 
         title = paste0("Polygenic Risk Score - ", title_suffix)) +
    scale_fill_manual(values = color_palette)
  
  # 保存图形
  mypath <- file.path("result", paste0("Boxplot PRS by NS group - ", title_suffix, ".png"))
  ggsave(mypath, plot = p, width = 10, height = 6, dpi = 300)
}


prs_boxplot(dta, "prs_multi_ancestry", "all NS groups", shpal)


















