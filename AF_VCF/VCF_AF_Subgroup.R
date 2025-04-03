
#subgroup
#Author: Haodi Yang 
#Date 7.6


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

# no.hom.0 = number of homozygotes for ref allele (00)
# no.het.01 = number of heterozygotes (01)
# no.het.10 = number of heterozygotes (10)
# no.het = number of heterozygotes (01 and 10)
# no.hom.1 = number of homozygotes for alt allele (11)

GTs$no.hom.0 <-  rowSums(GTs[,c(6:ncol(data))] == "0|0")   
GTs$no.het.01 <- rowSums(GTs[,c(6:ncol(data))] == "0|1")  
GTs$no.het.10 <- rowSums(GTs[,c(6:ncol(data))] == "1|0") 
GTs$no.het <- GTs$no.het.01+GTs$no.het.10
GTs$no.hom.1 <- rowSums(GTs[,c(6:ncol(data))] == "1|1") 
GTs$no.missing <- rowSums(GTs[,c(6:ncol(data))] == ".|.")  
GTs$total <- GTs$no.hom.0 + GTs$no.het + GTs$no.hom.1 + GTs$no.missing

GTs$alt.af <- (GTs$no.het + (GTs$no.hom.1*2))/(GTs$total*2)
GTs$alt.ac <- GTs$no.het + (GTs$no.hom.1*2)
GTs$an <- GTs$total*2

# Look at SNP statistics:
GTs[,c(1:5,(ncol(GTs)-9):ncol(GTs))]

###############chunk2: devide into 7 groups #######

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



###############chunk2: To subset by disease subgroups:

groupdata <- read.csv("data/Phenotype groups for exome data_SHclean.csv", header=T)

groupCFD <- subset(groupdata, groupdata$subgroup == "CFD")
groupnotry <- subset(groupdata, groupdata$subgroup == "INS - Steroids NOT Tried")
groupprisr <- subset(groupdata, groupdata$subgroup == "Primary Steroid Resistance")
groupparsr <- subset(groupdata, groupdata$subgroup == "Partial SR")
groupssns <- subset(groupdata, groupdata$subgroup == "SSNS")
groupsrnsgn <- subset(groupdata, groupdata$subgroup == "SRNS genetic")
groupsrnsnogn <- subset(groupdata, groupdata$subgroup == "SRNS non-genetic")


GTs.NoTry <- GTs[, c(1:5, which(names(GTs) %in% groupnotry$radar.ID))]
GTs.CFD <- GTs[, c(1:5, which(names(GTs) %in% groupCFD$radar.ID))]
GTs.PriSR <- GTs[, c(1:5, which(names(GTs) %in% groupprisr$radar.ID))]
GTs.ParSR <- GTs[,c(1:5, which(names(GTs) %in% groupparsr$radar.ID))]
GTs.SSNS <- GTs[,c(1:5, which(names(GTs) %in% groupssns$radar.ID))]
GTs.SRgene <- GTs[,c(1:5, which(names(GTs) %in% groupsrnsgn$radar.ID))]
GTs.SRnogn <- GTs[,c(1:5, which(names(GTs) %in% groupsrnsnogn$radar.ID))]


##################chunk2.1: CFD
GTs.CFD$no.hom.0 <-  rowSums(GTs.CFD[,c(1:ncol(GTs.CFD))] == "0|0")  
GTs.CFD$no.het.01 <- rowSums(GTs.CFD[,c(1:ncol(GTs.CFD))] == "0|1") 
GTs.CFD$no.het.10 <- rowSums(GTs.CFD[,c(1:ncol(GTs.CFD))] == "1|0")  
GTs.CFD$no.het <- GTs.CFD$no.het.01+GTs.CFD$no.het.10
GTs.CFD$no.hom.1 <- rowSums(GTs.CFD[,c(1:ncol(GTs.CFD))] == "1|1")  
GTs.CFD$no.missing <- rowSums(GTs.CFD[,c(1:ncol(GTs.CFD))] == ".|.")  
GTs.CFD$total <- GTs.CFD$no.hom.0 + GTs.CFD$no.het + GTs.CFD$no.hom.1 + GTs.CFD$no.missing

# alternative allele frequency (alt.af)
GTs.CFD$alt.af <- (GTs.CFD$no.het + (GTs.CFD$no.hom.1*2))/(GTs.CFD$total*2)

# alternative allele count (alt.ac)
GTs.CFD$alt.ac <- GTs.CFD$no.het + (GTs.CFD$no.hom.1*2)

# allele number (total number of patients * 2)
GTs.CFD$an <- GTs.CFD$total*2
GTs.CFD$PT <- "GTs.CFD"

CFD <- GTs.CFD[,c(1:5, (ncol(GTs.CFD)-10):ncol(GTs.CFD))]

write.table(CFD, "result/subgroup/CFD.txt", quote=F, row.names=F, col.names=T)
write.csv(CFD, "result/subgroup/CFD.csv", quote=F, row.names=F)

##################chunk2.2: NotTry
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
GTs.NoTry$PT <- "GTs.NoTry"

NotTry <- GTs.NoTry[,c(1:5, (ncol(GTs.NoTry)-10):ncol(GTs.NoTry))]

write.table(NotTry, "result/subgroup/NotTry.txt", quote=F, row.names=F, col.names=T)
write.csv(NotTry, "result/subgroup/NotTry.csv", quote=F, row.names=F)

##################chunk2.3: Primary SR
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
GTs.PriSR$PT <- "GTs.PriSR"

PrimarySR <- GTs.PriSR[,c(1:5, (ncol(GTs.PriSR)-10):ncol(GTs.PriSR))]

write.table(PrimarySR, "result/subgroup/PrimarySR.txt", quote=F, row.names=F, col.names=T)
write.csv(PrimarySR, "result/subgroup/PrimarySR.csv", quote=F, row.names=F)

##################chunk2.4: Partial SR
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
GTs.ParSR$PT <- "GTs.ParSR"

PartialSR <- GTs.ParSR[,c(1:5, (ncol(GTs.ParSR)-10):ncol(GTs.ParSR))]

write.table(PartialSR, "result/subgroup/PartialSR.txt", quote=F, row.names=F, col.names=T)
write.csv(PartialSR, "result/subgroup/PartialSR.csv", quote=F, row.names=F)

##################chunk2.5: SRNS Genetic
GTs.SRgene$no.hom.0 <-  rowSums(GTs.SRgene[,c(1:ncol(GTs.SRgene))] == "0|0")  
GTs.SRgene$no.het.01 <- rowSums(GTs.SRgene[,c(1:ncol(GTs.SRgene))] == "0|1") 
GTs.SRgene$no.het.10 <- rowSums(GTs.SRgene[,c(1:ncol(GTs.SRgene))] == "1|0")  
GTs.SRgene$no.het <- GTs.SRgene$no.het.01+GTs.SRgene$no.het.10
GTs.SRgene$no.hom.1 <- rowSums(GTs.SRgene[,c(1:ncol(GTs.SRgene))] == "1|1")  
GTs.SRgene$no.missing <- rowSums(GTs.SRgene[,c(1:ncol(GTs.SRgene))] == ".|.")  
GTs.SRgene$total <- GTs.SRgene$no.hom.0 + GTs.SRgene$no.het + GTs.SRgene$no.hom.1 + GTs.SRgene$no.missing

# alternative allele frequency (alt.af)
GTs.SRgene$alt.af <- (GTs.SRgene$no.het + (GTs.SRgene$no.hom.1*2))/(GTs.SRgene$total*2)

# alternative allele count (alt.ac)
GTs.SRgene$alt.ac <- GTs.SRgene$no.het + (GTs.SRgene$no.hom.1*2)

# allele number (total number of patients * 2)
GTs.SRgene$an <- GTs.SRgene$total*2
GTs.SRgene$PT <- "GTs.SRgene"

SRNS_Genetic <- GTs.SRgene[,c(1:5, (ncol(GTs.SRgene)-10):ncol(GTs.SRgene))]

write.table(SRNS_Genetic, "result/subgroup/SRNS_Genetic.txt", quote=F, row.names=F, col.names=T)
write.csv(SRNS_Genetic, "result/subgroup/SRNS_Genetic.csv", quote=F, row.names=F)

##################chunk2.6: SRNS non Genetic
GTs.SRnogn$no.hom.0 <-  rowSums(GTs.SRnogn[,c(1:ncol(GTs.SRnogn))] == "0|0")  
GTs.SRnogn$no.het.01 <- rowSums(GTs.SRnogn[,c(1:ncol(GTs.SRnogn))] == "0|1") 
GTs.SRnogn$no.het.10 <- rowSums(GTs.SRnogn[,c(1:ncol(GTs.SRnogn))] == "1|0")  
GTs.SRnogn$no.het <- GTs.SRnogn$no.het.01+GTs.SRnogn$no.het.10
GTs.SRnogn$no.hom.1 <- rowSums(GTs.SRnogn[,c(1:ncol(GTs.SRnogn))] == "1|1")  
GTs.SRnogn$no.missing <- rowSums(GTs.SRnogn[,c(1:ncol(GTs.SRnogn))] == ".|.")  
GTs.SRnogn$total <- GTs.SRnogn$no.hom.0 + GTs.SRnogn$no.het + GTs.SRnogn$no.hom.1 + GTs.SRnogn$no.missing

# alternative allele frequency (alt.af)
GTs.SRnogn$alt.af <- (GTs.SRnogn$no.het + (GTs.SRnogn$no.hom.1*2))/(GTs.SRnogn$total*2)

# alternative allele count (alt.ac)
GTs.SRnogn$alt.ac <- GTs.SRnogn$no.het + (GTs.SRnogn$no.hom.1*2)

# allele number (total number of patients * 2)
GTs.SRnogn$an <- GTs.SRnogn$total*2
GTs.SRnogn$PT <- "GTs.SRnogn"

SRNS_non_Genetic <- GTs.SRnogn[,c(1:5, (ncol(GTs.SRnogn)-10):ncol(GTs.SRnogn))]

write.table(SRNS_non_Genetic, "result/subgroup/SRNS_non_Genetic.txt", quote=F, row.names=F, col.names=T)
write.csv(SRNS_non_Genetic, "result/subgroup/SRNS_non_Genetic.csv", quote=F, row.names=F)

##################chunk2.7: SSNS
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

write.table(SSNS, "result/subgroup/SSNS.txt", quote=F, row.names=F, col.names=T)
write.csv(SSNS, "result/subgroup/SSNS.csv", quote=F, row.names=F)







