
library(dplyr)
library(readxl)
library(tidyverse)

rm(list=ls())


groupdata <- read.csv("data/Phenotype groups for exome data_SHclean.csv", header=T)

ethdata <- read.csv("data/Ethnicity groups for exome data.csv",header=T)

eth_group <- merge(groupdata, ethdata, by = "radar.ID", all = TRUE)

eth_group <- eth_group[eth_group$group != "Exclude patient", ]

group_eth_counts <- eth_group %>%
  group_by(group, nurture_ethnicity) %>%
  count() %>%
  arrange(group, nurture_ethnicity)

subgroup_eth_counts <- eth_group %>%
  group_by(subgroup, nurture_ethnicity) %>%
  count() %>%
  arrange(subgroup, nurture_ethnicity)


table(eth_group$nurture_ethnicity)


write.csv(group_eth_counts, "group_eth_counts.csv")
write.csv(subgroup_eth_counts, "subgroup_eth_counts.csv")

tableg <- group_eth_counts %>%
  pivot_wider(names_from = group, values_from = n, values_fill = list(n = 0))



tablesub <- subgroup_eth_counts %>%
  pivot_wider(names_from = subgroup, values_from = n, values_fill = list(n = 0))

write.csv(tableg, "table_group.csv")
write.csv(tablesub, "table_subgroup.csv")










