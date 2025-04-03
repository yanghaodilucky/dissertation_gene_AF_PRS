# Chi-square test
# Author: Haodi 

rm(list=ls())

library(readxl)
library(dplyr)

data <- read_excel("data/AF_result.xlsx")

######################data process########
data <- data %>%
  mutate(across(-SNPid, ~ as.integer(. * 1000)))

##################public vs whole group exon #########
public_all_exon <- data %>%
  rowwise() %>%
  mutate(
    Chi_Square = chisq.test(matrix(c(Public, 1000 - Public, Whole_group, 1000 - Whole_group), nrow = 2), correct = TRUE)$statistic,
    P_Value = chisq.test(matrix(c(Public, 1000 - Public, Whole_group, 1000 - Whole_group), nrow = 2), correct = TRUE)$p.value
  ) %>%
  ungroup() %>%  
  select(SNPid, Chi_Square, P_Value)

write.csv(public_all_exon, "result/chi/public_all_exon_results.csv")

##################public vs whole group not exon #########
public_noexon <- data %>%
  rowwise() %>%
  mutate(
    Chi_Square = chisq.test(matrix(c(Public, 1000 - Public, Whole_group_noexon, 1000 - Whole_group_noexon), nrow = 2), correct = TRUE)$statistic,
    P_Value = chisq.test(matrix(c(Public, 1000 - Public, Whole_group_noexon, 1000 - Whole_group_noexon), nrow = 2), correct = TRUE)$p.value
  ) %>%
  ungroup() %>%  
  select(SNPid, Chi_Square, P_Value)

write.csv(public_noexon, "result/chi/public_noexon_results.csv")





##################### public vs  group INS steroids not try####
public_gnottry <- data %>%
  rowwise() %>%
  mutate(
    Chi_Square = chisq.test(matrix(c(Public, 1000 - Public, g_nottry, 1000 - g_nottry), nrow = 2), correct = TRUE)$statistic,
    P_Value = chisq.test(matrix(c(Public, 1000 - Public, g_nottry, 1000 - g_nottry), nrow = 2), correct = TRUE)$p.value
  ) %>%
  ungroup() %>%  
  select(SNPid, Chi_Square, P_Value)

write.csv(public_gnottry, "result/chi/public_gnottry_results.csv")


################# public vs  group SSNS######
public_gssns <- data %>%
  rowwise() %>%
  mutate(
    Chi_Square = chisq.test(matrix(c(Public, 1000 - Public, g_ssns, 1000 - g_ssns), nrow = 2), correct = TRUE)$statistic,
    P_Value = chisq.test(matrix(c(Public, 1000 - Public, g_ssns, 1000 - g_ssns), nrow = 2), correct = TRUE)$p.value
  ) %>%
  ungroup() %>%  
  select(SNPid, Chi_Square, P_Value)

write.csv(public_gssns, "result/chi/public_gssns_results.csv")

################# public vs  group Primary SR######
public_gprisr <- data %>%
  rowwise() %>%
  mutate(
    Chi_Square = chisq.test(matrix(c(Public, 1000 - Public, g_prisr, 1000 - g_prisr), nrow = 2), correct = TRUE)$statistic,
    P_Value = chisq.test(matrix(c(Public, 1000 - Public, g_prisr, 1000 - g_prisr), nrow = 2), correct = TRUE)$p.value
  ) %>%
  ungroup() %>%  
  select(SNPid, Chi_Square, P_Value)

write.csv(public_gprisr, "result/chi/public_g_primarySR_results.csv")

################# public vs  group Secondary SR######
public_gsecsr <- data %>%
  rowwise() %>%
  mutate(
    Chi_Square = chisq.test(matrix(c(Public, 1000 - Public, g_secsr, 1000 - g_secsr), nrow = 2), correct = TRUE)$statistic,
    P_Value = chisq.test(matrix(c(Public, 1000 - Public, g_secsr, 1000 - g_secsr), nrow = 2), correct = TRUE)$p.value
  ) %>%
  ungroup() %>%  
  select(SNPid, Chi_Square, P_Value)

write.csv(public_gsecsr, "result/chi/public_g_secondarySR_results.csv")

################# public vs  group Secondary SR######
public_gparsr <- data %>%
  rowwise() %>%
  mutate(
    Chi_Square = chisq.test(matrix(c(Public, 1000 - Public, g_parsr, 1000 - g_parsr), nrow = 2), correct = TRUE)$statistic,
    P_Value = chisq.test(matrix(c(Public, 1000 - Public, g_parsr, 1000 - g_parsr), nrow = 2), correct = TRUE)$p.value
  ) %>%
  ungroup() %>%  
  select(SNPid, Chi_Square, P_Value)

write.csv(public_gparsr, "result/chi/public_g_partialSR_results.csv")

################# public vs  subgroup CFD######
public_subCFD <- data %>%
  rowwise() %>%
  mutate(
    Chi_Square = chisq.test(matrix(c(Public, 1000 - Public, sub_cfd, 1000 - sub_cfd), nrow = 2), correct = TRUE)$statistic,
    P_Value = chisq.test(matrix(c(Public, 1000 - Public, sub_cfd, 1000 - sub_cfd), nrow = 2), correct = TRUE)$p.value
  ) %>%
  ungroup() %>%  
  select(SNPid, Chi_Square, P_Value)

write.csv(public_subCFD, "result/chi/public_subgroup_CFD_results.csv")

################# public vs  subgroup not try######
public_subnottry <- data %>%
  rowwise() %>%
  mutate(
    Chi_Square = chisq.test(matrix(c(Public, 1000 - Public, sub_nottry, 1000 - sub_nottry), nrow = 2), correct = TRUE)$statistic,
    P_Value = chisq.test(matrix(c(Public, 1000 - Public, sub_nottry, 1000 - sub_nottry), nrow = 2), correct = TRUE)$p.value
  ) %>%
  ungroup() %>%  
  select(SNPid, Chi_Square, P_Value)

write.csv(public_subnottry, "result/chi/public_subgroup_nottry_results.csv")


################# public vs  subgroup primary SR######
public_subprisr <- data %>%
  rowwise() %>%
  mutate(
    Chi_Square = chisq.test(matrix(c(Public, 1000 - Public, sub_prisr, 1000 - sub_prisr), nrow = 2), correct = TRUE)$statistic,
    P_Value = chisq.test(matrix(c(Public, 1000 - Public, sub_prisr, 1000 - sub_prisr), nrow = 2), correct = TRUE)$p.value
  ) %>%
  ungroup() %>%  
  select(SNPid, Chi_Square, P_Value)

write.csv(public_subprisr, "result/chi/public_subgroup_primarySR_results.csv")

################# public vs  subgroup partial SR######
public_subparsr <- data %>%
  rowwise() %>%
  mutate(
    Chi_Square = chisq.test(matrix(c(Public, 1000 - Public, sub_parsr, 1000 - sub_parsr), nrow = 2), correct = TRUE)$statistic,
    P_Value = chisq.test(matrix(c(Public, 1000 - Public, sub_parsr, 1000 - sub_parsr), nrow = 2), correct = TRUE)$p.value
  ) %>%
  ungroup() %>%  
  select(SNPid, Chi_Square, P_Value)

write.csv(public_subparsr, "result/chi/public_subgroup_partialSR_results.csv")

################# public vs  subgroup SR gene######
public_subsrgene <- data %>%
  rowwise() %>%
  mutate(
    Chi_Square = chisq.test(matrix(c(Public, 1000 - Public, sub_srge, 1000 - sub_srge), nrow = 2), correct = TRUE)$statistic,
    P_Value = chisq.test(matrix(c(Public, 1000 - Public, sub_srge, 1000 - sub_srge), nrow = 2), correct = TRUE)$p.value
  ) %>%
  ungroup() %>%  
  select(SNPid, Chi_Square, P_Value)

write.csv(public_subsrgene, "result/chi/public_subgroup_SR_gene_results.csv")


################# public vs  subgroup SR non gene######
public_subsrnoge <- data %>%
  rowwise() %>%
  mutate(
    Chi_Square = chisq.test(matrix(c(Public, 1000 - Public, sub_srnoge, 1000 - sub_srnoge), nrow = 2), correct = TRUE)$statistic,
    P_Value = chisq.test(matrix(c(Public, 1000 - Public, sub_srnoge, 1000 - sub_srnoge), nrow = 2), correct = TRUE)$p.value
  ) %>%
  ungroup() %>%  
  select(SNPid, Chi_Square, P_Value)

write.csv(public_subsrnoge, "result/chi/public_subgroup_SR_non_gene_results.csv")

################# public vs  subgroup SSNS######
public_subssns <- data %>%
  rowwise() %>%
  mutate(
    Chi_Square = chisq.test(matrix(c(Public, 1000 - Public, sub_ssns, 1000 - sub_ssns), nrow = 2), correct = TRUE)$statistic,
    P_Value = chisq.test(matrix(c(Public, 1000 - Public, sub_ssns, 1000 - sub_ssns), nrow = 2), correct = TRUE)$p.value
  ) %>%
  ungroup() %>%  
  select(SNPid, Chi_Square, P_Value)

write.csv(public_subssns, "result/chi/public_subgroup_SSNS_results.csv")



