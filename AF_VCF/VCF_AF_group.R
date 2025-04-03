# This script loads in SNP data (VCF format) and calculates the SNP Allele frequency for each SNP in the file (allele frequency = allele count/allele number)

library(dplyr)
library(readxl)

rm(list=ls())

##############chunk 1: all patients###################
data <- read.table("data/all_chrs.GWAS_ET_NS_2023_Sam144NS_Rasheed_INS_HaodiRiskSNPs_GTs_R2_0.9.vcf",comment.char = "",skip=42,header=T)

# To remove repeating elements of participant 'SentrixPosition'/Sample ID column names (e.g. the first few repeating characters)-so that they then match with phenotype files later:
newnames1 <- gsub("^.{0,34}", "", names(data[,c(10:473)]))
newnames2 <- gsub("^.{0,42}", "", names(data[,c(474:482)]))
newnames3 <- gsub("^.{0,46}", "", names(data[,c(483:570)]))
newnames4 <- gsub("^.{0,50}", "", names(data[,c(571:614)]))

newnames5 <- ""
for (i in 615:ncol(data)){
	newnames5[i-614] <- gsub("B", "", unlist(strsplit(names(data[,i,drop=F]), "_"))[1])
}

newnames <- c(newnames1, newnames2, newnames3, newnames4, newnames5)

data <- data[,c(1:5,10:ncol(data))]
names(data) <- c("chr","pos","RSID","ref","alt",newnames)
GTs <- data



##########chunk1.1: change colnames to radar.ID##################

# first I change sampleIDs into RADAR id to group them
radarid1 <- read.csv("data/GSA_sample_sheet_Amy_forCCA.csv", header=T)

radarid2 <- read_excel("data/NS_GSA_samples_with_withdrawn_detailAO.xlsx")

radarid2 <- radarid2[-(1:6), ]

radarid2 <- radarid2 %>% rename(radar_ID = PPID)

# deal with the radarid1

ra1 <- radarid1[, c("Sample_ID", "radar_ID")]
vecra1 <- setNames(as.character(ra1$radar_ID), as.character(ra1$Sample_ID))
new_colnames <- sapply(names(GTs), function(colname) {
  if (colname %in% names(vecra1)) {
    return(vecra1[[colname]])
  } else {
    return(colname)
  }
})

names(GTs) <- new_colnames



# deal with radarid2
ra2 <- radarid2[, c("SentrixPosition", "radar_ID")]
vecra2 <- setNames(as.character(ra2$radar_ID), as.character(ra2$SentrixPosition))
new_colnames <- sapply(names(GTs), function(colname) {
  if (colname %in% names(vecra2)) {
    return(vecra2[[colname]])
  } else {
    return(colname)
  }
})

names(GTs) <- new_colnames


##########chunk1.2: calculate AF for exon or not exon patients####
groupdata <- read.csv("data/Phenotype groups for exome data_SHclean.csv", header=T)

ethdata <- read.csv("data/Ethnicity groups for exome data.csv",header=T)

GTs.exon <- GTs[, c(1:5, which(names(GTs) %in% groupdata$radar.ID))]

GTs_columns <- names(GTs)
cols_to_include <- setdiff(GTs_columns, names(GTs.exon))
cols_to_include <- c(GTs_columns[1:5], cols_to_include)
cols_to_include <- na.omit(cols_to_include)
GTs.noexon <- GTs[, cols_to_include]


# no.hom.0 = number of homozygotes for ref allele (00)
# no.het.01 = number of heterozygotes (01)
# no.het.10 = number of heterozygotes (10)
# no.het = number of heterozygotes (01 and 10)
# no.hom.1 = number of homozygotes for alt allele (11)


# calculate variants within exon location
GTs.exon$no.hom.0 <-  rowSums(GTs.exon[,c(6:ncol(GTs.exon))] == "0|0")   
GTs.exon$no.het.01 <- rowSums(GTs.exon[,c(6:ncol(GTs.exon))] == "0|1")  
GTs.exon$no.het.10 <- rowSums(GTs.exon[,c(6:ncol(GTs.exon))] == "1|0") 
GTs.exon$no.het <- GTs.exon$no.het.01+GTs.exon$no.het.10
GTs.exon$no.hom.1 <- rowSums(GTs.exon[,c(6:ncol(GTs.exon))] == "1|1") 
GTs.exon$no.missing <- rowSums(GTs.exon[,c(6:ncol(GTs.exon))] == ".|.")  
GTs.exon$total <- GTs.exon$no.hom.0 + GTs.exon$no.het + GTs.exon$no.hom.1 + GTs.exon$no.missing

GTs.exon$alt.af <- (GTs.exon$no.het + (GTs.exon$no.hom.1*2))/(GTs.exon$total*2)
GTs.exon$alt.ac <- GTs.exon$no.het + (GTs.exon$no.hom.1*2)
GTs.exon$an <- GTs.exon$total*2


all_AF_exon <- GTs.exon[,c(1:5, (ncol(GTs.exon)-9):ncol(GTs.exon))]

write.table(all_AF_exon, "result/group/all_patient_exon.txt", row.names = F)
write.csv(all_AF_exon, "result/group/all_patient_exon.csv", row.names = F)


# calculate variants not in exon locations

GTs.noexon$no.hom.0 <-  rowSums(GTs.noexon[,c(6:ncol(GTs.noexon))] == "0|0")   
GTs.noexon$no.het.01 <- rowSums(GTs.noexon[,c(6:ncol(GTs.noexon))] == "0|1")  
GTs.noexon$no.het.10 <- rowSums(GTs.noexon[,c(6:ncol(GTs.noexon))] == "1|0") 
GTs.noexon$no.het <- GTs.noexon$no.het.01+GTs.noexon$no.het.10
GTs.noexon$no.hom.1 <- rowSums(GTs.noexon[,c(6:ncol(GTs.noexon))] == "1|1") 
GTs.noexon$no.missing <- rowSums(GTs.noexon[,c(6:ncol(GTs.noexon))] == ".|.")  
GTs.noexon$total <- GTs.noexon$no.hom.0 + GTs.noexon$no.het + GTs.noexon$no.hom.1 + GTs.noexon$no.missing

GTs.noexon$alt.af <- (GTs.noexon$no.het + (GTs.noexon$no.hom.1*2))/(GTs.noexon$total*2)
GTs.noexon$alt.ac <- GTs.noexon$no.het + (GTs.noexon$no.hom.1*2)
GTs.noexon$an <- GTs.noexon$total*2


all_AF_noexon <- GTs.noexon[,c(1:5, (ncol(GTs.noexon)-9):ncol(GTs.noexon))]

write.table(all_AF_noexon, "result/group/all_patient_noexon.txt", row.names = F)
write.csv(all_AF_noexon, "result/group/all_patient_noexon.csv", row.names = F)




###############chunk2.0: To subset by disease groups:##############

groupnotry <- subset(groupdata, groupdata$group == "INS - Steroids NOT Tried")
groupprisr <- subset(groupdata, groupdata$group == "Primary Steroid Resistance")
groupsecsr <- subset(groupdata, groupdata$group == "Secondary Steroid Resistance")
groupparsr <- subset(groupdata, groupdata$group == "Partial Steroid Resistance")
groupssns <- subset(groupdata, groupdata$group == "SSNS")

GTs.NoTry <- GTs[, c(1:5, which(names(GTs) %in% groupnotry$radar.ID))]
GTs.PriSR <- GTs[, c(1:5, which(names(GTs) %in% groupprisr$radar.ID))]
GTs.SecSR <- GTs[,c(1:5, which(names(GTs) %in% groupsecsr$radar.ID))]
GTs.ParSR <- GTs[,c(1:5, which(names(GTs) %in% groupparsr$radar.ID))]
GTs.SSNS <- GTs[,c(1:5, which(names(GTs) %in% groupssns$radar.ID))]


###########chunk2.1: group1 SSNS###########
GTs.SSNS$no.hom.0 <-  rowSums(GTs.SSNS[,c(1:ncol(GTs.SSNS))] == "0|0")  
GTs.SSNS$no.het.01 <- rowSums(GTs.SSNS[,c(1:ncol(GTs.SSNS))] == "0|1") 
GTs.SSNS$no.het.10 <- rowSums(GTs.SSNS[,c(1:ncol(GTs.SSNS))] == "1|0")  
GTs.SSNS$no.het <- GTs.SSNS$no.het.01+GTs.SSNS$no.het.10
GTs.SSNS$no.hom.1 <- rowSums(GTs.SSNS[,c(1:ncol(GTs.SSNS))] == "1|1")  
GTs.SSNS$no.missing <- rowSums(GTs.SSNS[,c(1:ncol(GTs.SSNS))] == ".|.")  
GTs.SSNS$total <- GTs.SSNS$no.hom.0 + GTs.SSNS$no.het + GTs.SSNS$no.hom.1 + GTs.SSNS$no.missing

# alternative allele frequency (alt.af)
GTs.SSNS$alt.af <- (GTs.SSNS$no.het + (GTs.SSNS$no.hom.1*2))/(GTs.SSNS$total*2)

# alternative allele count (alt.ac)
GTs.SSNS$alt.ac <- GTs.SSNS$no.het + (GTs.SSNS$no.hom.1*2)

# allele number (total number of patients * 2)
GTs.SSNS$an <- GTs.SSNS$total*2
GTs.SSNS$PT <- "GTs.SSNS"

SSNS <- GTs.SSNS[,c(1:5, (ncol(GTs.SSNS)-10):ncol(GTs.SSNS))]

write.table(SSNS, "result/group/SSNS.txt", quote=F, row.names=F, col.names=T)
write.csv(SSNS, "result/group/SSNS.csv", quote=F, row.names=F)


###########chunk2.2: group2:no try###########
GTs.NoTry$no.hom.0 <-  rowSums(GTs.NoTry[,c(1:ncol(GTs.NoTry))] == "0|0")  
GTs.NoTry$no.het.01 <- rowSums(GTs.NoTry[,c(1:ncol(GTs.NoTry))] == "0|1") 
GTs.NoTry$no.het.10 <- rowSums(GTs.NoTry[,c(1:ncol(GTs.NoTry))] == "1|0")  
GTs.NoTry$no.het <- GTs.NoTry$no.het.01+GTs.NoTry$no.het.10
GTs.NoTry$no.hom.1 <- rowSums(GTs.NoTry[,c(1:ncol(GTs.NoTry))] == "1|1")  
GTs.NoTry$no.missing <- rowSums(GTs.NoTry[,c(1:ncol(GTs.NoTry))] == ".|.")  
GTs.NoTry$total <- GTs.NoTry$no.hom.0 + GTs.NoTry$no.het + GTs.NoTry$no.hom.1 + GTs.NoTry$no.missing

# alternative allele frequency (alt.af)
GTs.NoTry$alt.af <- (GTs.NoTry$no.het + (GTs.NoTry$no.hom.1*2))/(GTs.NoTry$total*2)

# alternative allele count (alt.ac)
GTs.NoTry$alt.ac <- GTs.NoTry$no.het + (GTs.NoTry$no.hom.1*2)

# allele number (total number of patients * 2)
GTs.NoTry$an <- GTs.NoTry$total*2
GTs.NoTry$PT <- "INS_Sterios not try"

INS_Steroids_not_try <- GTs.NoTry[,c(1:5, (ncol(GTs.NoTry)-10):ncol(GTs.NoTry))]

write.table(INS_Steroids_not_try, "result/group/INS_Steroids_not_try.txt", quote=F, row.names=F, col.names=T)
write.csv(INS_Steroids_not_try, "result/group/INS_Steroids_not_try.csv", quote=F, row.names=F)



#############chunk2.3: group3: Primary Steroid Resistance ###############

GTs.PriSR$no.hom.0 <-  rowSums(GTs.PriSR[,c(1:ncol(GTs.PriSR))] == "0|0")  
GTs.PriSR$no.het.01 <- rowSums(GTs.PriSR[,c(1:ncol(GTs.PriSR))] == "0|1") 
GTs.PriSR$no.het.10 <- rowSums(GTs.PriSR[,c(1:ncol(GTs.PriSR))] == "1|0")  
GTs.PriSR$no.het <- GTs.PriSR$no.het.01+GTs.PriSR$no.het.10
GTs.PriSR$no.hom.1 <- rowSums(GTs.PriSR[,c(1:ncol(GTs.PriSR))] == "1|1")  
GTs.PriSR$no.missing <- rowSums(GTs.PriSR[,c(1:ncol(GTs.PriSR))] == ".|.")  
GTs.PriSR$total <- GTs.PriSR$no.hom.0 + GTs.PriSR$no.het + GTs.PriSR$no.hom.1 + GTs.PriSR$no.missing

# alternative allele frequency (alt.af)
GTs.PriSR$alt.af <- (GTs.PriSR$no.het + (GTs.PriSR$no.hom.1*2))/(GTs.PriSR$total*2)

# alternative allele count (alt.ac)
GTs.PriSR$alt.ac <- GTs.PriSR$no.het + (GTs.PriSR$no.hom.1*2)

# allele number (total number of patients * 2)
GTs.PriSR$an <- GTs.PriSR$total*2
GTs.PriSR$PT <- "PriSR"

Primary_steroid_resistance <- GTs.PriSR[,c(1:5, (ncol(GTs.PriSR)-10):ncol(GTs.PriSR))]

write.table(Primary_steroid_resistance, "result/group/Primary_steroid_resistance.txt", quote=F, row.names=F, col.names=T)
write.csv(Primary_steroid_resistance, "result/group/Primary_steroid_resistance.csv", quote=F, row.names=F)

###############chunk2.4: group4:Secondary Steroid Resistance##############

GTs.SecSR$no.hom.0 <-  rowSums(GTs.SecSR[,c(1:ncol(GTs.SecSR))] == "0|0")  
GTs.SecSR$no.het.01 <- rowSums(GTs.SecSR[,c(1:ncol(GTs.SecSR))] == "0|1") 
GTs.SecSR$no.het.10 <- rowSums(GTs.SecSR[,c(1:ncol(GTs.SecSR))] == "1|0")  
GTs.SecSR$no.het <- GTs.SecSR$no.het.01+GTs.SecSR$no.het.10
GTs.SecSR$no.hom.1 <- rowSums(GTs.SecSR[,c(1:ncol(GTs.SecSR))] == "1|1")  
GTs.SecSR$no.missing <- rowSums(GTs.SecSR[,c(1:ncol(GTs.SecSR))] == ".|.")  
GTs.SecSR$total <- GTs.SecSR$no.hom.0 + GTs.SecSR$no.het + GTs.SecSR$no.hom.1 + GTs.SecSR$no.missing

# alternative allele frequency (alt.af)
GTs.SecSR$alt.af <- (GTs.SecSR$no.het + (GTs.SecSR$no.hom.1*2))/(GTs.SecSR$total*2)

# alternative allele count (alt.ac)
GTs.SecSR$alt.ac <- GTs.SecSR$no.het + (GTs.SecSR$no.hom.1*2)

# allele number (total number of patients * 2)
GTs.SecSR$an <- GTs.SecSR$total*2
GTs.SecSR$PT <- "SecSR"

Secondery_steroid_resistance <- GTs.SecSR[,c(1:5, (ncol(GTs.SecSR)-10):ncol(GTs.SecSR))]

write.table(Secondery_steroid_resistance, "result/group/Secondery_steroid_resistance.txt", quote=F, row.names=F, col.names=T)
write.csv(Secondery_steroid_resistance, "result/group/Secondery_steroid_resistance.csv", quote=F, row.names=F)

##################chunk2:5 group5: Partial Steroid Resistance

GTs.ParSR$no.hom.0 <-  rowSums(GTs.ParSR[,c(1:ncol(GTs.ParSR))] == "0|0")  
GTs.ParSR$no.het.01 <- rowSums(GTs.ParSR[,c(1:ncol(GTs.ParSR))] == "0|1") 
GTs.ParSR$no.het.10 <- rowSums(GTs.ParSR[,c(1:ncol(GTs.ParSR))] == "1|0")  
GTs.ParSR$no.het <- GTs.ParSR$no.het.01+GTs.ParSR$no.het.10
GTs.ParSR$no.hom.1 <- rowSums(GTs.ParSR[,c(1:ncol(GTs.ParSR))] == "1|1")  
GTs.ParSR$no.missing <- rowSums(GTs.ParSR[,c(1:ncol(GTs.ParSR))] == ".|.")  
GTs.ParSR$total <- GTs.ParSR$no.hom.0 + GTs.ParSR$no.het + GTs.ParSR$no.hom.1 + GTs.ParSR$no.missing

# alternative allele frequency (alt.af)
GTs.ParSR$alt.af <- (GTs.ParSR$no.het + (GTs.ParSR$no.hom.1*2))/(GTs.ParSR$total*2)

# alternative allele count (alt.ac)
GTs.ParSR$alt.ac <- GTs.ParSR$no.het + (GTs.ParSR$no.hom.1*2)

# allele number (total number of patients * 2)
GTs.ParSR$an <- GTs.ParSR$total*2
GTs.ParSR$PT <- "ParSR"

Partial_steroid_resistance <- GTs.ParSR[,c(1:5, (ncol(GTs.ParSR)-10):ncol(GTs.ParSR))]

write.table(Partial_steroid_resistance, "result/group/Partial_steroid_resistance.txt", quote=F, row.names=F, col.names=T)
write.csv(Partial_steroid_resistance, "result/group/Partial_steroid_resistance.csv", quote=F, row.names=F)



