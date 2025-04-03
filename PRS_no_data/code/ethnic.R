rm(list=ls())


library(ggplot2)
library(dplyr)
library(tidyr)

ethnic_data <- data.frame(
  Group = c("Asian (Chinese)",
            "Asian (Pakistani, Indian or Bangladeshi)",
            "Black ",
            "Mixed ethnicity",
            "White"),
  Min = c(1.50387028, -0.98870955, -0.16086967, -0.32617013, -1.08050809),
  Q1 = c(1.60251569, -0.42998093, 0.56441700, 0.01394916, -0.19053219),
  Median = c(1.70116109, 0.03118129, 0.90582904, 0.47247199, 0.21823884),
  Mean = c(1.70116109, 0.02621116, 0.91033954, 0.58913146, 0.22009505),
  Q3 = c(1.79980650, 0.43665343, 1.21671369, 0.95365527, 0.63511868),
  Max = c(1.89845190, 1.01185820, 1.91047582, 1.91055126, 1.68756896)
)

# 创建箱线图
boxplot <- ggplot(ethnic_data, aes(x = Group, fill = Group)) +
  geom_boxplot(aes(ymin = Min, lower = Q1, middle = Median, upper = Q3, ymax = Max),
               stat = "identity") +
  scale_fill_manual(values = c("#FF9999", "#66B2FF", "#99FF99", "#FFCC99", "#FF9966")) + # 自定义颜色
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10), 
    axis.text.y = element_text(face = "bold"), 
    axis.text.x = element_text(face = "bold", angle = 0, hjust = 0.5, vjust = 0.5),  # x轴标签居中
    axis.title.x = element_blank(), 
    axis.title.y = element_text(face = "bold"), 
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(), # 隐藏纵向的网格线
    panel.grid.minor.x = element_blank(), # 隐藏次纵向的网格线
    panel.grid.major.y = element_line(color = "grey80"), # 添加水平网格线
    panel.grid.minor.y = element_blank(), # 隐藏次水平网格线
    axis.line = element_line(color = "black") # 显示轴线
  ) +
  labs(title = "Boxplot of PrS in Different Ethnic Groups", 
       x = "Ethnic Group", 
       y = "PRS") +
  coord_cartesian(clip = "off")  # 保持竖直方向

ggsave(filename = "result/PRS_vs_ethnic.png", plot = boxplot, width = 13, height = 8, dpi = 300)
