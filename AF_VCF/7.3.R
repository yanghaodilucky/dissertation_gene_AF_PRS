
# auhor: Haodi
# data 2024.7.2
# to calculate AF from VCF document


library(dplyr)

rm(list=ls())

##########chunk 1: read the data and clean it and change the colume names##########
dta <- read.table(
  "NURTUREINS_Sam144NS_RasheedINS_TopMed_AllAnc_AllChr_Barry2023_SNPs_R2_0.9_ACs_NoMissingData.vcf",
  comment.char = "",skip=44,header=T
  )

# To remove repeating elements of participant 'SentrixPosition'/Sample ID column names 
# (e.g. the first few repeating characters)-so that they then match with phenotype files later:
newnames1 <- gsub("^.{0,42}", "", names(dta[,c(10:18)]))
newnames2 <- gsub("^.{0,46}", "", names(dta[,c(19:106)]))
newnames3 <- gsub("^.{0,50}", "", names(dta[,c(107:150)]))
newnames4 <- gsub("^.{0,34}", "", names(dta[,c(151:205)]))
newnames5 <- gsub("B(\\d+)_.*", "\\1",names(dta[,c(206:269)]))  

newnames <- c(newnames1, newnames2, newnames3, newnames4, newnames5)

dta <- dta[,c(1:5,10:ncol(dta))]  
names(dta) <- c("chr","pos","RSID","ref","alt",newnames)

rsid <- c("rs62397901", "rs2746432", "rs73885319", "rs60910145", "rs71785313",
          "rs487575", "rs59882675", "rs2076523", "rs9348883", "rs2637678",
          "rs2858829", "rs76615866", "rs190705792", "rs55730955", "rs8062322",
          "rs4649032", "rs12431424", "rs10817678", "rs2858829", "rs111796602",
          "rs1801274", "rs16946160", "rs9303279", "rs530462", "rs2844580", 
          "rs755622", "rs751611955", "rs12911841", "rs2857607", "rs2285450",
          "rs412175", "rs6531527", "rs62397901", "rs115180879", "rs10518133", 
          "rs1805732", "rs113752715", "rs1264705", "rs12431424", "rs6584128", 
          "rs4649032", "rs11086243", "rs4979462", "rs10817678", "rs7848647", 
          "rs1012507", "rs1265889", "rs751611955", "rs2637678", "rs1264705", 
          "rs6916716", "rs9261947", "rs9348883", "rs2076523", "rs6763024", 
          "rs2674382", "rs1799937", "rs2301254", "rs6508", "rs115180879", 
          "rs34217742", "rs1891621")
uniqrsid <- unique(rsid)
uniqrsid


dta <- subset(dta, RSID%in%uniqrsid)



GTs <- dta


#############chunk2: calculate all patients' AF#####################
# no.hom.0 = number of homozygotes for ref allele (0)
# no.het.01 = number of heterozygotes (01)
# no.het.10 = number of heterozygotes (10)
# no.het = number of heterozygotes (01 and 10) (1)
# no.hom.1 = number of homozygotes for alt allele (2)

GTs$no.hom.0 <-  rowSums(GTs[,c(6:ncol(dta))] == "0")   
GTs$no.het <- rowSums(GTs[,c(6:ncol(dta))] == "1")
GTs$no.hom.1 <- rowSums(GTs[,c(6:ncol(dta))] == "2") 
GTs$no.missing <- rowSums(GTs[,c(6:ncol(dta))] == ".")  
GTs$total <- GTs$no.hom.0 + GTs$no.het + GTs$no.hom.1 + GTs$no.missing

GTs$alt.af <- (GTs$no.het + (GTs$no.hom.1*2))/(GTs$total*2)
GTs$alt.ac <- GTs$no.het + (GTs$no.hom.1*2)
GTs$an <- GTs$total*2

# Look at SNP statistics:
GTs[,c(1:5,(ncol(GTs)-9):ncol(GTs))]

# output as txt and csv format
all_AF <- GTs[,c(1:5, (ncol(GTs)-9):ncol(GTs))]

write.table(all_AF, "all_patient_AF.txt", row.names = F)
write.csv(all_AF, "all_patient_AF.csv", row.names = F)





