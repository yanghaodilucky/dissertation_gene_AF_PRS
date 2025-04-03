## Purpose: Compare PRSs between the different NS subgroups
## Input: NS_DNAm_cohort_SSNS_PRS_Barry.csv and NS_clin_lab_med_phenotypes.csv
## Output: 

# Notes: Only 260 patients from the DNAm NS cohort have SNP data
# Only 38% of the PRS SNPs are present in the NS SNP data


#-----------------------------------------------------------------
# Load packages and set up directories
#-----------------------------------------------------------------        
my_packages <- c("devtools", "readxl", "EpiFunc", "janitor", "lubridate", "skimr", 
                 "gtsummary", "ggplot2", "naniar", "openxlsx", 
                 "multcomp", "tidyverse", "survival") 
lapply(my_packages, require, character.only = TRUE)

# Set directories
project.dir <- "~/1. Work/My work/Epigenetics - INS"
data.dir <- file.path(project.dir, "Analysis", "DNAm PhD data analysis")
results.dir <- file.path(project.dir, "Results", "Thesis data", "PRS")

#-----------------------------------------------------------------
# Load data and combine clinical and prs data
#-----------------------------------------------------------------        
prs <- read.csv(file.path(data.dir, "NS_DNAm_cohort_SSNS_PRS_Barry.csv"))
prs$patient_id <- as.character(prs$patient_id)

clin <- read.csv(file.path(data.dir, "NS_clin_lab_med_phenotypes.csv")) # clinical subgroup data
clin$rough_NSgroup <- as.factor(clin$rough_NSgroup)

clin_prs <- left_join(prs, clin, by = "patient_id") # 75 patients missing PRS information
rm(clin, prs)

# subset the european patients to be used for the euro-prs
euro_pts_only <- clin_prs[clin_prs$ethnicity_group == "White", ]

# Select just the 4 subgroups of interest, all ancestries
grps <- c("SSNS", "SRNS genetic", "SRNS non-genetic", "CFD")
four_nsgrps <- clin_prs[clin_prs$rough_NSgroup %in% grps,]

# Select just the 4 subgroups of interest and the European patients
euro_nsgrps <- clin_prs[clin_prs$rough_NSgroup %in% grps & clin_prs$ethnicity_group == "White",]

#-----------------------------------------------------------------
# Compare PRS scores across the 4 NS subgroups
#-----------------------------------------------------------------        

results <- aggregate(clin_prs$prs_multi_ancestry, list(clin_prs$rough_NSgroup), FUN=summary) 
#write.csv(results, file.path(results.dir, "PRS by NS subgroup summary.csv"), row.names = FALSE)

# Anova to test for differences between the groups
a_res <- aov(prs_multi_ancestry ~ rough_NSgroup, data = four_nsgrps) 
summary(a_res)
# Which exact groups differ?
# Tukey HSD test:
post_test <- glht(a_res,linfct = mcp(rough_NSgroup = "Tukey"))
summary(post_test)

# Just euro patients with euro PRS
a_res_euro <- aov(prs_euro_ancestry ~ rough_NSgroup, data = euro_nsgrps) 
summary(a_res_euro)
post_test_euro <- glht(a_res_euro,linfct = mcp(rough_NSgroup = "Tukey"))
summary(post_test_euro)

#-----------------------------------------------------------------
# Plots
#-----------------------------------------------------------------        

# Box plot - PRS values by each subgroup
shpal <- c("deeppink4", "rosybrown3", "midnightblue", "lightblue4", "darkorange2", "darkturquoise", "aquamarine4", "wheat3", "darkslategray", "midnightblue")
shpal2 <- c("deeppink4", "midnightblue", "lightblue4","darkorange2")

prs_plot <- function(df1, y1, z1, sh1) {
  p <- df1 %>% ggplot(aes(x = reorder(rough_NSgroup, y1, decreasing=T), y = y1, fill = rough_NSgroup)) + 
    geom_boxplot(width=0.5) +
    theme(plot.title=element_text(hjust=0.5, size = 10), axis.text.y = element_text(face="bold"), 
          axis.text.x = element_text(face="bold", angle = -70 ), axis.title.x = element_blank(), 
          axis.title.y = element_text(face="bold"), legend.title=element_blank()) +
    labs(y = "Polygenic risk score ", title = paste0("SSNS polygenic risk score - ", z1)) +
    scale_fill_manual(values = sh1)
  mypath <- file.path(results.dir, paste0("Boxplot SSNS PRS by NS group - ", z1, ".png"))
  png(file = mypath) 
  plot(p)
  dev.off() 
}

# All NS subgroups
prs_plot(clin_prs, clin_prs$prs_multi_ancestry, "multi ancestry PRS - all pts", shpal)
#prs_plot(clin_prs, clin_prs$prs_euro_ancestry, "euro ancestry PRS - all pts", shpal)
prs_plot(euro_pts_only, euro_pts_only$prs_euro_ancestry, "euro ancestry PRS - euro pts", shpal)

# 4 NS subgroups only 
prs_plot(four_nsgrps, four_nsgrps$prs_multi_ancestry, "multi ancestry PRS - 4 subgroups", shpal2)
prs_plot(euro_nsgrps, euro_nsgrps$prs_euro_ancestry, "euro ancestry PRS - euro pts, 4 subgroups", shpal2)

# scatterplot - PRS values by each subgroup
prs_scatter <- function(df1, y1, z1, sh1) {
  p <- df1 %>% ggplot(aes(x = reorder(rough_NSgroup, y1, decreasing=T), y = y1, fill = rough_NSgroup)) + 
    geom_jitter(aes(colour = rough_NSgroup), 
                width = 0.15, size = 3) +
    theme(plot.title=element_text(hjust=0.5, size = 10), axis.text.y = element_text(face="bold"), 
          axis.text.x = element_text(face="bold", angle = -70 ), axis.title.x = element_blank(), 
          axis.title.y = element_text(face="bold"), legend.title=element_blank()) +
    labs(y = "Polygenic risk score ", title = paste0("SSNS polygenic risk score - ", z1)) +
    scale_colour_manual(values = sh1)
  mypath <- file.path(results.dir, paste0("Scatterplot SSNS PRS by NS group - ", z1, ".png"))
  png(file = mypath) 
  plot(p)
  dev.off() 
}

# All NS subgroups
prs_scatter(clin_prs, clin_prs$prs_multi_ancestry, "multi ancestry PRS - all pts", shpal)
prs_scatter(euro_pts_only, euro_pts_only$prs_euro_ancestry, "euro ancestry PRS - euro pts", shpal)

# 4 NS subgroups only 
prs_scatter(four_nsgrps, four_nsgrps$prs_multi_ancestry, "multi ancestry PRS - 4 subgroups", shpal2)
prs_scatter(euro_nsgrps, euro_nsgrps$prs_euro_ancestry, "euro ancestry PRS - euro pts, 4 subgroups", shpal2)
