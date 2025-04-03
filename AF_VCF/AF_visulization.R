# visulization
# author: Haodi Yang
# date: 12.7.2024

# this script is the visulization of the AF result 


library(dplyr)
library(readxl)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)


rm(list=ls())

data <- read_excel("data/AF_result.xlsx")

set1_colors <- brewer.pal(3, "Set1")
colors <- c(set1_colors[3], set1_colors[1], set1_colors[2]) # green red blue



########## the first plot: public and all patients(exon and noexon) ########


data1 <- melt(data, id.vars = "SNPid", 
                  measure.vars = c("Public", "Whole_group", "Whole_group_noexon"),
                  variable.name = "Group", value.name = "AF")


public_patients <- ggplot(data1, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS groups",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


  
ggsave("visulization/Public_whole_group.pdf", plot = public_patients, width = 24, height = 8)


############plot2: public and only exon variants###########


data2 <- melt(data, id.vars = "SNPid", 
                  measure.vars = c("Public", "Whole_group"),
                  variable.name = "Group", value.name = "AF")


public_exon <- ggplot(data2, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS groups (only exon variants)",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggsave("visulization/Public_exon.pdf", plot = public_exon, width = 24, height = 8)



############plot3: public and other variants from SNP data###########


data3 <- melt(data, id.vars = "SNPid", 
              measure.vars = c("Public", "Whole_group_noexon"),
              variable.name = "Group", value.name = "AF")


public_notexon <- ggplot(data3, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS groups (other variants from SNP data)",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggsave("visulization/Public_not_exon.pdf", plot = public_notexon, width = 24, height = 8)


############plot4: public and group: INS-Steroids not try###########


data4 <- melt(data, id.vars = "SNPid", 
              measure.vars = c("Public", "g_nottry"),
              variable.name = "Group", value.name = "AF")


public_g_nottry <- ggplot(data4, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS group: INS-Steroids not try",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggsave("visulization/Public_g_nottry.pdf", plot = public_g_nottry, width = 24, height = 8)


############plot5: public and group: SSNS###########


data5 <- melt(data, id.vars = "SNPid", 
              measure.vars = c("Public", "g_ssns"),
              variable.name = "Group", value.name = "AF")


public_g_ssns <- ggplot(data5, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS group: SSNS",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggsave("visulization/Public_g_ssns.pdf", plot = public_g_ssns, width = 24, height = 8)

############plot6: public and group: Primary SR###########


data6 <- melt(data, id.vars = "SNPid", 
              measure.vars = c("Public", "g_prisr"),
              variable.name = "Group", value.name = "AF")


public_g_prisr <- ggplot(data6, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS group: Primary SR",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggsave("visulization/Public_g_prisr.pdf", plot = public_g_prisr, width = 24, height = 8)



############plot7: public and group: Secondary SR###########


data7 <- melt(data, id.vars = "SNPid", 
              measure.vars = c("Public", "g_secsr"),
              variable.name = "Group", value.name = "AF")


public_g_secsr <- ggplot(data7, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS group: Secondary SR",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggsave("visulization/Public_g_secsr.pdf", plot = public_g_secsr, width = 24, height = 8)

############plot8: public and group: Partial SR###########


data8 <- melt(data, id.vars = "SNPid", 
              measure.vars = c("Public", "g_parsr"),
              variable.name = "Group", value.name = "AF")


public_g_parsr <- ggplot(data8, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS group: Partial SR",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggsave("visulization/Public_g_parsr.pdf", plot = public_g_parsr, width = 24, height = 8)

############plot9: public and subgroup: CFD###########


data9 <- melt(data, id.vars = "SNPid", 
              measure.vars = c("Public", "sub_cfd"),
              variable.name = "Group", value.name = "AF")


public_sub_cfd <- ggplot(data9, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS subgroup: CFD",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggsave("visulization/Public_sub_cfd.pdf", plot = public_sub_cfd, width = 24, height = 8)

############plot10: public and subgroup: primary SR###########


data10 <- melt(data, id.vars = "SNPid", 
              measure.vars = c("Public", "sub_cfd"),
              variable.name = "Group", value.name = "AF")


public_sub_prisr <- ggplot(data10, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS subgroup: Primary SR",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggsave("visulization/Public_sub_prisr.pdf", plot = public_sub_prisr, width = 24, height = 8)


############plot11: public and subgroup: Partial SR###########


data11 <- melt(data, id.vars = "SNPid", 
               measure.vars = c("Public", "sub_parsr"),
               variable.name = "Group", value.name = "AF")


public_sub_parsr <- ggplot(data11, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS subgroup: Partial SR",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggsave("visulization/Public_sub_parsr.pdf", plot = public_sub_parsr, width = 24, height = 8)

############plot12: public and subgroup: SRNS Genetic###########


data12 <- melt(data, id.vars = "SNPid", 
               measure.vars = c("Public", "sub_srge"),
               variable.name = "Group", value.name = "AF")


public_sub_srge <- ggplot(data12, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS subgroup: SRNS Genetic",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggsave("visulization/Public_sub_srge.pdf", plot = public_sub_srge, width = 24, height = 8)


############plot13: public and subgroup: SRNS non-Genetic###########


data13 <- melt(data, id.vars = "SNPid", 
               measure.vars = c("Public", "sub_srnoge"),
               variable.name = "Group", value.name = "AF")


public_sub_srnoge <- ggplot(data13, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS subgroup: SRNS non-Genetic",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggsave("visulization/Public_sub_srnoge.pdf", plot = public_sub_srnoge, width = 24, height = 8)

############plot14: public and subgroup: SSNS###########


data14 <- melt(data, id.vars = "SNPid", 
               measure.vars = c("Public", "sub_ssns"),
               variable.name = "Group", value.name = "AF")


public_sub_ssns <- ggplot(data14, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS subgroup: SSNS",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggsave("visulization/Public_sub_ssns.pdf", plot = public_sub_ssns, width = 24, height = 8)

############plot15: public and subgroup: INS-Steroids not try###########


data15 <- melt(data, id.vars = "SNPid", 
               measure.vars = c("Public", "sub_nottry"),
               variable.name = "Group", value.name = "AF")


public_sub_nottry <- ggplot(data15, aes(x = SNPid, y = AF, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = colors) +
  labs(title = "AF Values Comparison between public and NS subgroup: INS-Steroids not try",
       x = "SNPID",
       y = "AF Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)) 
theme(axis.text.x = element_text(angle = 90, hjust = 1))



ggsave("visulization/Public_sub_nottry.pdf", plot = public_sub_nottry, width = 24, height = 8)

######################## combine the plots together ######################

comb_public_exon_nonexon <- grid.arrange(public_patients, public_exon, public_notexon, ncol = 1)


ggsave("visulization/combine_exon_nonexon.pdf", plot = comb_public_exon_nonexon, width = 24, height = 24)



com_group <- grid.arrange(public_g_nottry, public_g_ssns, public_g_prisr, 
                           public_g_secsr, public_g_parsr,
                           ncol = 1)


ggsave("visulization/combine_groups.pdf", plot = com_group, width = 30, height = 40)

com_subgroup <- grid.arrange(public_sub_cfd, public_sub_prisr, public_sub_parsr, 
                             public_sub_srge, public_sub_srnoge, public_sub_ssns, 
                             public_sub_nottry, 
                          ncol = 1)


ggsave("visulization/combine_subgroups.pdf", plot = com_subgroup, width = 30, height = 49)





