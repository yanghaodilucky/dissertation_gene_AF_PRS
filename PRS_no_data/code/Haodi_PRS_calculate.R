# calculate PRS for 260 patients
# Haodi 
## Input: Allele counts for all NS patients
## Output: PRS for all NS patients

rm(list=ls())
ns_snps <- read.table("data/prs/NURTUREINS_Sam144NS_RasheedINS_TopMed_AllAnc_AllChr_Barry2023_SNPs_R2_0.9_ACs_NoMissingData.vcf",comment.char = "",skip=44,header=T)


my_packages <- c( "readxl", "EpiFunc", "janitor", "lubridate", 
                  "skimr", "gtsummary", "ggplot2", "naniar", "readr",
                  "openxlsx", "tidyverse", "survival", "utils", "Cairo") 

lapply(my_packages, require, character.only = TRUE)

multi_prs <- read.table("data/prs/Multipop_PRS_TopmedID_hg38_EURtrained.txt.gz", header = TRUE) # 1,081,086 SNPS

# Clean PRS data so that it matches NS SNP data
clean_ns_snp <- function(df){
  df$CHROM <- gsub(pattern = ":.*", replacement="", df$TopmedID)
  df$CHROM <- gsub(pattern = "chr", replacement="", df$CHROM)
  df$POS <- gsub(pattern = "chr[0-9]:", replacement="", df$TopmedID)
  df$POS <- gsub(pattern = ":.*", replacement="", df$POS)
  df$link <- paste0(df$CHROM, "_", df$POS)
  return(df)
}


multi_prs <- clean_ns_snp(multi_prs)

# CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT,
names(ns_snps)[names(ns_snps) == "X.CHROM"] <- "CHROM"
ns_snps$link <- paste0(ns_snps$CHROM, "_", ns_snps$POS)

# Check which PRS SNPs are in the NS dataset
available_snps <- ns_snps$link
multi_prs$available <- if_else(multi_prs$link %in% available_snps, "Yes", "No")


################# change patients name to radar ID #####################

newnames1 <- sub(".*(SH[0-9]{8}).*", "\\1", names(ns_snps[,c(10:150)]))
newnames2 <- sub(".*_(\\d+_R\\d+C\\d+)_.*", "\\1", names(ns_snps[,c(151:205)]))
newnames3 <- sub("B(\\d+)_.*", "\\1", names(ns_snps[,c(206:269)]))

newnames <- c(newnames1, newnames2, newnames3)

ns_snps <- ns_snps[,c(1:5,10:ncol(ns_snps))]
names(ns_snps) <- c("CHROM","POS","ID","REF","ALT",newnames, "link")

radarid1 <- read.csv("data/prs/GSA_sample_sheet_Amy_forCCA.csv", header=T)
radarid2 <- read_excel("data/prs/NS_GSA_samples_with_withdrawn_detailAO.xlsx")
radarid2 <- radarid2[-(1:6), ]
radarid2 <- radarid2 %>% rename(radar_ID = PPID)

gts_colnames <- colnames(ns_snps)

for (i in seq_along(gts_colnames)) {
  sample_id <- gts_colnames[i]
  
  radar_id <- radarid1$radar_ID[radarid1$Sample_ID == sample_id]
  
  
  if (length(radar_id) > 0) {
    gts_colnames[i] <- as.character(radar_id)
  }
}

unique_colnames <- make.unique(gts_colnames)
colnames(ns_snps) <- unique_colnames


gts_colnames <- colnames(ns_snps)
for (i in seq_along(gts_colnames)) {
  sentrixposition <- gts_colnames[i]
  radar_id <- radarid2$radar_ID[radarid2$SentrixPosition == sentrixposition]
  if (length(radar_id) > 0) {
    gts_colnames[i] <- as.character(radar_id)
  }
}

unique_colnames <- make.unique(gts_colnames)

colnames(ns_snps) <- unique_colnames

################## calculate the PRS #################


m_prs <- multi_prs[multi_prs$available == "Yes", ]
multi_ns <- left_join(m_prs, ns_snps, by = "link")

multi_weight <- multi_ns[, c("TopmedID", "weight")]
rownames(multi_weight) <- multi_weight$TopmedID
rownames(multi_ns) <- multi_ns$TopmedID

# Select only numeric columns for SNP data
genotype_cols <- setdiff(colnames(ns_snps), c("CHROM", "POS", "ID", "REF", "ALT", "link"))
snp_data <- multi_ns[, genotype_cols]
snp_data <- snp_data %>% mutate(across(everything(), as.numeric))

# Match the order of the coefficients to the order of the SNPs
coefs <- as.numeric(multi_weight$weight)
score <- coefs * as.matrix(snp_data)
prs_multi_ancestry <- colSums(score)

prs_multi <- as.data.frame(prs_multi_ancestry)
prs_multi$name <- as.character(rownames(prs_multi))

write.csv(prs_multi,"result/prs.csv", row.names = TRUE)












