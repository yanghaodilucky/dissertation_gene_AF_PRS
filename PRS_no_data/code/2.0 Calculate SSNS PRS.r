## Purpose: Calculate the SSNS polygenic risk score (Barry et al) for the NS patients in the DNAm cohort
## Input: Allele counts for all NS patients
## Output: PRS for all NS patients

# Notes: Only 260 patients from the DNAm NS cohort have SNP data
# Only 38% of the PRS SNPs are present in the NS SNP data


#-----------------------------------------------------------------
# Load packages and set up directories
#-----------------------------------------------------------------
my_packages <- c( "readxl", "EpiFunc", "janitor", "lubridate", 
                 "skimr", "gtsummary", "ggplot2", "naniar", "readr",
                 "openxlsx", "tidyverse", "survival", "utils", "Cairo") 

lapply(my_packages, require, character.only = TRUE)

# Set directories

project.dir <- "/projects/MRC-IEU/research/data/radar/epigenetic/ins/methylation" # if run on epi-franklin
data.dir <- file.path(project.dir, "dev", "release_candidate", "data", "output-20240318")
raw.dir <- file.path(project.dir, "raw")
snp.dir <- file.path(raw.dir, "NS polygenic risk score SNPs")

output.dir.data <- file.path(data.dir, "ML Elastic net models")
scripts.dir <- file.path(project.dir, "dev", "release_candidate", "scripts", "mlr3")

#-----------------------------------------------------------------
# Unzip files
#-----------------------------------------------------------------
# Unzip and load file (uses utils package)
# Specify the path to the Zip file
zip_file <- file.path(snp.dir,"SSNS_PRS_Barry2023.zip")

# Extract files from the zip archive and save in the snp directory
unzip(zip_file, exdir = snp.dir)
# Files are:
# Multipop_PRS_TopmedID_hg38_EURtrained.txt.gz
# Euro_only_PRS_TopmedID_hg38.txt.gz

#-----------------------------------------------------------------
# Load PRS and NS SNP data
#-----------------------------------------------------------------

# Load PRSs (Barry et al 2023)
multi_prs <- read.table(file.path(snp.dir,"Multipop_PRS_TopmedID_hg38_EURtrained.txt.gz"), header = TRUE) # 1,081,086 SNPS
euro_prs <- read.table(file.path(snp.dir,"Euro_only_PRS_TopmedID_hg38.txt.gz"), header = TRUE) # 1,075,893 SNPS

# Number of SNPs in common - 1,075,893
both <- intersect(multi_prs$TopmedID, euro_prs$TopmedID)
length(both) 

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
euro_prs <- clean_ns_snp(euro_prs)
# TopmedID; A1; weight; CHROM; POS


# Load the NS data
ns_snps <- read.table(file.path(snp.dir,"NURTUREINS_Sam144NS_RasheedINS_TopMed_AllAnc_AllChr_Barry2023_SNPs_R2_0.9_ACs_NoMissingData.vcf.gz"), header=TRUE, skip=44, comment.char = "" )
names(ns_snps)[names(ns_snps) == "X.CHROM"] <- "CHROM"
# CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT,
# cols names are exome ids but with X in front of any numbers
ns_snps$link <- paste0(ns_snps$CHROM, "_", ns_snps$POS)

# Check which PRS SNPs are in the NS dataset
available_snps <- ns_snps$link
multi_prs$available <- if_else(multi_prs$link %in% available_snps, "Yes", "No")
euro_prs$available <- if_else(euro_prs$link %in% available_snps, "Yes", "No")

# Graph of available and missing SNPS with their coefficients
#shpal2 <- c("darkorange2", "darkslategray")
#p.width <- 1500
#p.height <- 500

#CairoPNG(file.path(output.dir.data, "SSNS Barry et al PRS available SNPS.png"),
#         width=p.width, height=p.height)
#p <- ggplot(data = multi_prs, aes(y = weight, x = reorder(TopmedID, weight, decreasing=T), color = available)) +
#  geom_point(size = 0.5) + 
  # scale_color_brewer(palette = "dark2") +
#  scale_color_manual(values = shpal2) +
#  geom_hline(yintercept=0, color='dark red', linetype='dashed', alpha=.5) +
#  labs(title='SNPs in the PRS
#       ', 
#       x='
#       ', y = 'PRS weighting
#       ') +
#  theme(plot.title=element_text(hjust=0.5, size = 14, face = "bold"), 
#        axis.text.x = element_blank(), axis.text.y = element_text(size = 12), 
#        axis.title.x = element_text(size = 12, face="bold"), axis.title.y = element_text(size = 12, face="bold"),
#        legend.position = "right", 
#        panel.background = element_rect(color="dark grey", fill="grey95"),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank())
#print(p)
#dev.off() 

# Just select available SNPs and join the PRS and NS SNPs:
m_prs <- multi_prs[multi_prs$available == "Yes",]
multi_ns <- left_join(m_prs, ns_snps, by = "link")

e_prs <- euro_prs[euro_prs$available == "Yes",]
euro_ns <- left_join(e_prs, ns_snps, by = "link")

#-----------------------------------------------------------------
# Load linker file for Radar IDs and NS SNP data = 263 NS patients have SNP data that passed QC
#-----------------------------------------------------------------

# Load the linker IDs - SNP:Radar
ids <- read_excel(file.path(snp.dir, "NS_330_IDs_genomic_barcodes_v3.xlsx"))
ids <- as.data.frame(ids)
names(ids)[names(ids) == "sample name for extracting from imputed VCF"] <- "name"

# Remove samples without SNP data
ids <- ids[!is.na(ids$GSA_ID),]

# Remove samples with dodgy QC
fail_qc <- c("nephros_UK_Rasheed, but filtered out due to missing GT data", "Rasheed - removed at QC due to missing GT rate >10%")
ids <- ids[!ids$cohort %in% fail_qc,]
# radar id, GWAS_id, patient_id, cohort, name

# Add X to the names of any generated SNP data from IEU so that it matches the name of the SNP files
ids$name <- if_else(ids$cohort == "Sam_144", paste0("X",ids$name), ids$name)
# name in the ID file links to the column names in multi_ns and euro_ns

#-----------------------------------------------------------------
# Calculate multi-ancestry PRS from SNP data
#-----------------------------------------------------------------
# Df with rownames as snps and column of weightings
multi_weight <- multi_ns[, c("TopmedID", "weight")] # 411,512 SNPS
rownames(multi_weight) <- multi_weight$TopmedID
# 13 duplicate SNPs are in there but with different alternative alleles - not an issue as patients will just have allele counts of 0 for the ones they don't have. 

# Make a matrix with rownames as snps and columns as patients (with X)
snp_ids <- ids$name # just take the patients that the SNP data passed the QC
rownames(multi_ns) <- multi_ns$TopmedID
multi_ns <- multi_ns[, colnames(multi_ns) %in% snp_ids]

# Calculate multi-ancestry PRS
# Check only taking sites from the PRS (not necessary but keep in in case I use this code for other reasons)
sites <- intersect(rownames(multi_ns), rownames(multi_weight))
multi_ns <- multi_ns[sites,,drop=F] # produce a vector of SNP counts (,, is a matrix thing)

# Match the order of the coefficients to the order of the SNPs
multi_weight <- multi_weight[match(rownames(multi_ns), rownames(multi_weight)),] 

# Calculate PRS
coefs <- multi_weight$weight
score <- coefs * multi_ns
prs_multi_ancestry <- colSums(score) # column names are patient ids

# Add the PRS to the clinical data
prs_multi <- as.data.frame(prs_multi_ancestry)
prs_multi$name <- as.character(rownames(prs_multi))
clin_prs <- left_join(ids, prs_multi, by = "name")
clin_multi_prs <- clin_prs[, c("patient_id", "prs_multi_ancestry", "name")]

#-----------------------------------------------------------------
# Calculate European ancestry PRS from SNP data
#-----------------------------------------------------------------
euro_weight <- euro_ns[, c("TopmedID", "weight")] # 409,812 SNPS
rownames(euro_weight) <- euro_weight$TopmedID

rownames(euro_ns) <- euro_ns$TopmedID
euro_ns <- euro_ns[, colnames(euro_ns) %in% snp_ids]

# Calculate European-ancestry PRS
# Check only taking sites from the PRS
sites <- intersect(rownames(euro_ns), rownames(euro_weight))
euro_ns <- euro_ns[sites,,drop=F] # produce a vector of SNP counts (,, is a matrix thing)

# Match the order of the coefficients to the order of the SNPs
euro_weight <- euro_weight[match(rownames(euro_ns), rownames(euro_weight)),] 

# Calculate PRS
coefs <- euro_weight$weight
score <- coefs * euro_ns
prs_euro_ancestry <- colSums(score) # column names are patient ids

# Add the PRS to the clinical data
prs_euro <- as.data.frame(prs_euro_ancestry)
prs_euro$name <- as.character(rownames(prs_euro))
clin_prs_both <- left_join(clin_multi_prs, prs_euro, by = "name")
clin_prs_both <- clin_prs_both[, c("patient_id", "prs_multi_ancestry", "prs_euro_ancestry", "name")]

#-----------------------------------------------------------------
# Save - 260 patients with usable PRSs
#-----------------------------------------------------------------
# Remove any patients missing PRS
clin_prs_both <- clin_prs_both[!is.na(clin_prs_both$prs_multi_ancestry),]

# Save
write.csv(clin_prs_both, file = file.path(output.dir.data, "NS_DNAm_cohort_SSNS_PRS_Barry.csv"), row.names = FALSE)

#-----------------------------------------------------------------
# 3 patients without PRSs - Radar IDs: 351, 453, 525
#-----------------------------------------------------------------

missing_euro <- euro_ns[complete.cases(euro_ns),]

# Don't appear to have SNP data - investigate further if needed
check <- euro_ns[, c("X118_SH00000705_118_SH00000705_118_SH00000705_118_SH00000705",
                     "X77_SH00000840_77_SH00000840_77_SH00000840_77_SH00000840",
                     "X93_SH00000768_93_SH00000768_93_SH00000768_93_SH00000768")] 

check <- ns_snps[, c("X1_SH00000807_1_SH00000807_1_SH00000807_1_SH00000807")] # present in data

# No SNPs in the data from Amy
check <- ns_snps[, c("X118_SH00000705_118_SH00000705_118_SH00000705_118_SH00000705")] 
check <- ns_snps[, c("X77_SH00000840_77_SH00000840_77_SH00000840_77_SH00000840")] 
check <- ns_snps[, c("X93_SH00000768_93_SH00000768_93_SH00000768_93_SH00000768")] 


