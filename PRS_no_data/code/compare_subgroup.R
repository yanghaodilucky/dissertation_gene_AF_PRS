# subgroup PRS and Visulization

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


# compare between subgroups

dta <- dta[dta$rough_NSgroup != "Unclear gen/non - exclude in subgroup analyses",]
dta <- dta[!is.na(dta$first_diagnosis_age), ]
dta$rough_NSgroup <- as.factor(dta$rough_NSgroup)

grps <- c("SSNS", "SRNS genetic", "SRNS non-genetic", "CFD")
four_nsgrps <- dta[dta$rough_NSgroup %in% grps,]

results <- aggregate(dta$prs_multi_ancestry, list(dta$rough_NSgroup), FUN=summary)
write.csv(results, "result/PRS by NS subgroup summary.csv", row.names = FALSE)

a_res <- aov(prs_multi_ancestry ~ rough_NSgroup, data = dta) 
a_NS_res <- summary(a_res)

# Which exact groups differ?
# Tukey HSD test:
post_test <- glht(a_res,linfct = mcp(rough_NSgroup = "Tukey"))
post_NS_res <- summary(post_test)


library(dplyr)

dta <- dta %>%
  mutate(rough_NSgroup = recode(rough_NSgroup, 
                          "SRNS genetic" = "Monogenic SRNS",
                          "SRNS non-genetic" = "Non Monogenic SRNS"))


#-----------------------------------------------------------------
# Plots
#-----------------------------------------------------------------        


shpal <- c("navyblue", "orange", "green", "red", "purple", "brown", "cyan")
shpal2 <- c("navyblue", "green", "red", "purple")



prs_boxplot <- function(df, y_var, title_suffix, color_palette) {
  p <- df %>% 
    ggplot(aes(x = reorder(rough_NSgroup, !!sym(y_var), decreasing = TRUE), 
               y = !!sym(y_var), fill = rough_NSgroup)) + 
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
    labs(y = "Polygenic risk score", 
         title = paste0("SSNS polygenic risk score - ", title_suffix)) +
    scale_fill_manual(values = color_palette)+
   mypath <- file.path("result", paste0("Boxplot SSNS PRS by NS group - ", title_suffix, ".png"))
  
  
  ggsave(mypath, plot = p, width = 8, height = 6, dpi = 300)
}
prs_boxplot(dta, "prs_multi_ancestry", "all", shpal)
prs_boxplot(four_nsgrps, "prs_multi_ancestry", "4 subgroups", shpal2) 
  

prs_scatter <- function(df, y_var, title_suffix, color_palette) {
  p <- df %>% 
    ggplot(aes(x = reorder(rough_NSgroup, !!sym(y_var), decreasing = TRUE), 
               y = !!sym(y_var), colour = rough_NSgroup)) + 
    geom_jitter(width = 0.15, size = 3) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10), 
      axis.text.y = element_text(face = "bold"), 
      axis.text.x = element_text(face = "bold", angle = 0), 
      axis.title.x = element_blank(), 
      axis.title.y = element_text(face = "bold"), 
      legend.title = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major.x = element_blank(), # 隐藏纵向的网格线
      panel.grid.minor.x = element_blank(), # 隐藏次纵向的网格线
      panel.grid.major.y = element_line(color = "grey80"), # 添加水平网格线
      panel.grid.minor.y = element_blank(), # 隐藏次水平网格线
      axis.line = element_line(color = "black") # 显示轴线
    ) +
    labs(y = "Polygenic risk score", 
         title = paste0("SSNS polygenic risk score - ", title_suffix)) +
    scale_colour_manual(values = color_palette)
  
  
  mypath <- file.path("result", paste0("Scatterplot SSNS PRS by NS group - ", title_suffix, ".png"))
  
 
  ggsave(mypath, plot = p, width = 8, height = 6, dpi = 300)
}

prs_scatter(dta, "prs_multi_ancestry", "all", shpal)
prs_scatter(four_nsgrps, "prs_multi_ancestry", "4 subgroups", shpal2)
  

################ linear regression ##################

model1 <- lm(first_diagnosis_age ~ prs_multi_ancestry, data = dta)

summary(model1)
confint(model1)

cbind(model1$coefficients,confint(model1))

dta$linear_fitted <- predict(model1, dta)

plot1 <- ggplot(dta) +
  geom_point(mapping=aes(x=prs_multi_ancestry,y=first_diagnosis_age)) +
  geom_line(mapping=aes(x=prs_multi_ancestry,y=linear_fitted), color="blue",size=1.5)+
  labs(title = "Scatterplot of PRS against age with Fitted Values from Linear Regression", 
       x = "PRS", y = "age") 



# 
ggsave(filename = "result/Scatterplot_PRS_vs_Age.png", plot = plot1, width = 8, height = 6)


################ different models ##################
poly_model <- lm(first_diagnosis_age ~ poly(prs_multi_ancestry,3), data = dta)
summary(poly_model)
confint(poly_model)

cbind(poly_model$coefficients,confint(poly_model))

dta$poly_fitted <- predict(poly_model, dta)

library(splines)
spline_model <- lm(first_diagnosis_age ~ bs(prs_multi_ancestry, df = 5), data = dta)

summary(spline_model)
confint(spline_model)

cbind(spline_model$coefficients,confint(spline_model))

dta$spline_fitted <- predict(spline_model, dta)

library(mgcv)
gam_model <- gam(first_diagnosis_age ~ s(prs_multi_ancestry), data = dta)
summary(gam_model)
confint(gam_model)

cbind(gam_model$coefficients,confint(gam_model))

dta$gam_fitted <- predict(gam_model, dta)

plot2 <- ggplot(dta, aes(x = prs_multi_ancestry, y = first_diagnosis_age)) +
  geom_point() +
  geom_line(aes(y = linear_fitted, color = "Linear"), size = 1, linetype = "solid") +
  geom_line(aes(y = poly_fitted, color = "Polynomial"), size = 1, linetype = "solid") +
  geom_line(aes(y = spline_fitted, color = "Spline"), size = 1, linetype = "solid") +
  geom_line(aes(y = gam_fitted, color = "GAM"), size = 1, linetype = "solid") +
  labs(title = "Age vs PRS: Different Relationships",
       x = "PRS",
       y = "Age",
       color = "Model") +
  scale_color_manual(values = c("Linear" = "blue", "Polynomial" = "red", 
                                "Spline" = "pink", "GAM" = "cyan")) +
  theme_bw() +
  theme(legend.position = "right")


ggsave(filename = "result/Scatterplot_PRS_vs_Age_disserent_relationship.png", plot = plot2, width = 10, height = 8, dpi = 300)

####################### compare between different genders ##################
results_gender <- aggregate(dta$prs_multi_ancestry, list(dta$gender), FUN=summary)
write.csv(results_gender, "result/PRS by gender summary.csv", row.names = FALSE)

gender_res <- aov(prs_multi_ancestry ~ gender, data = dta) 
anova_gender_res <- summary(gender_res)


####################### compare between different ethical groups##################
dta_eth <- dta[dta$ethnicity_group != "Unknown",]
dta_eth$ethnicity_group <- as.factor(dta_eth$ethnicity_group)

results_eth <- aggregate(dta_eth$prs_multi_ancestry, list(dta_eth$ethnicity_group), FUN=summary)
write.csv(results_eth, "result/PRS by ethnicity summary.csv", row.names = FALSE)

eth_res <- aov(prs_multi_ancestry ~ ethnicity_group, data = dta_eth) 
anova_eth_res <- summary(eth_res)

post_eth_test <- glht(eth_res,linfct = mcp(ethnicity_group = "Tukey"))
post_eth_res <- summary(post_eth_test)


# output results

anova_eth_results <- capture.output(anova_eth_res)
writeLines(anova_eth_results, "result/ANOVA_ethnicity_results.txt")
tukey_eth_results <- capture.output(post_eth_res)
writeLines(tukey_eth_results, "result/Tukey_HSD_ethnicity_results.txt")

anova_gender_results <- capture.output(anova_gender_res)
writeLines(anova_gender_results, "result/ANOVA_gender_results.txt")



anova_group_results <- capture.output(a_NS_res)
writeLines(anova_group_results, "result/ANOVA_subgroup_results.txt")
tukey_group_results <- capture.output(post_NS_res)
writeLines(tukey_group_results, "result/Tukey_HSD_subgroup_results.txt")













