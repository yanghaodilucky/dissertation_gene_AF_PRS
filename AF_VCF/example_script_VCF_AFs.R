# This script loads in SNP data (VCF format) and calculates the SNP Allele frequency for each SNP in the file (allele frequency = allele count/allele number)

data <- read.table("~/forHaodi/all_chrs.GWAS_ET_NS_2023_Sam144NS_Rasheed_INS_HaodiRiskSNPs_GTs_R2_0.9.vcf",comment.char = "",skip=42,header=T)

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


# To subset by a disease group (or 'PT') - one example:
CKD_PTs <- read.csv("~/Sam_144_NS/GSA_sample_sheet_Amy_forCCA.csv", header=T)

subgroup1 <- CKD_PTs[CKD_PTs$Phenotype.Group == "SSNS",]
GTs.SSNS <- GTs[,names(GTs) %in% subgroup1$Sample_ID]

# no.hom.0 = number of homozygotes for ref allele (00)
# no.het.01 = number of heterozygotes (01)
# no.het.10 = number of heterozygotes (10)
# no.het = number of heterozygotes (01 and 10)
# no.hom.1 = number of homozygotes for alt allele (11)

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

write.table(GTs.SSNS, "filename.txt", quote=F, row.names=F, col.names=T)
