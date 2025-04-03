# backup

sh_ids <- c("SH00000816", "SH00000811", "SH00000814", "SH00000826", "SH00000818", "SH00000821", "SH00000824", "SH00000827",
            "SH00000870", "SH00000828", "SH00000865", "SH00000867", "SH00000753", "SH00000715", "SH00000734", "SH00000750",
            "SH00000755", "SH00000784", "SH00000831", "SH00000832", "SH00000726", "SH00000774", "SH00000834", "SH00000835",
            "SH00000712", "SH00000739", "SH00000838", "SH00000839", "SH00000787", "SH00000719", "SH00000729", "SH00000741",
            "SH00000777", "SH00000730", "SH00000742", "SH00000778", "SH00000844", "SH00000731", "SH00000769", "SH00000862",
            "SH00000722", "SH00000732", "SH00000738", "SH00000780", "SH00000850", "SH00000708", "SH00000772", "SH00000725",
            "SH00000736")

radar_ids <- c(2, 6, 9, 12, 13, 16, 19, 21, 22, 23, 26, 28, 32, 33, 37, 38, 42, 46, 47, 48, 50, 54, 57, 58, 59, 61, 68, 75, 76, 79,
               80, 81, 84, 90, 91, 94, 97, 100, 103, 105, 109, 110, 111, 114, 127, 128, 133, 139, 140)


mapping <- setNames(radar_ids, sh_ids)

original_colnames <- names(GTs)

new_colnames <- sapply(original_colnames, function(colname) {
  if (colname %in% names(mapping)) {
    return(as.character(mapping[[colname]]))
  } else {
    return(colname)
  }
})

names(GTs) <- new_colnames

write.csv(GTs, "with_radarid.csv")












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

#end of chunk 1

all_AF <- GTs[,c(1:5, (ncol(GTs)-9):ncol(GTs))]

write.table(all_AF, "result/group/all_patient_AF.txt", row.names = F)
write.csv(all_AF, "result/group/all_patient_AF.csv", row.names = F)
# To subset by a disease group: