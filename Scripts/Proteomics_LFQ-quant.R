# This script was created by Leandro Xavier Neves, PhD (https://orcid.org/0000-0002-6074-1025)

# Load packages
library("dplyr")
library("ggplot2")
library("reshape")
library("hrbrthemes")
library("viridis")

# Set working directory
setwd("Y:/Leandro/2024/Projects/18_T7capping")

# Load Fragpipe output files
data <- read.table("./01_timsTOF-HT_transfection_1/Fragpipe/LN_LFQ-MBR_trypsinKRnotP_1MC_skyline_noMBR/combined_protein.tsv",
                   sep = "\t", header = TRUE)

# Subset MNV proteins
MNV <- dplyr::filter(data, Organism == "Norovirus (isolate Mouse/")

# Subset and log2 transform quantitative data
quant_columns <- colnames(MNV[,210:274])
quant_columns <- sort(quant_columns)

int <- MNV[,c("Description",quant_columns)]
int[,quant_columns][int[,quant_columns] == 0] <- NA
int <- dplyr::mutate(int,across(where(is.numeric),log2))

# Plot box-plot
conditions <- c(rep(paste0("Cond_0", 2:9), each = 5),rep(paste0("Cond_", 10:13), each = 5),rep("Mock", each = 5))
conditions <- conditions[c(61:65,1:60)]

int_long <- melt(int)
int_long <- int_long[c(61:65,1:60),]
int_long <- cbind(int_long,conditions)
int_long$conditions <- factor(int_long$conditions, levels = unique(conditions)) #converts the conditions column to a factor and specifies the desired order of the levels

Fig2B <- ggplot(int_long, aes(x = conditions, y = value, fill = conditions)) +
  #geom_boxplot() +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.2) +
  #geom_boxplot(coef=NULL) + 
  stat_boxplot(geom = "errorbar", width = 0.3, coef = NULL, alpha = 0.8, linewidth = 0.2) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  scale_y_continuous(limits = c(16,22)) +
  #geom_jitter(color="black", size=0.6, alpha=0.9, width = 0, shape = 1, stroke = 0.2) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0) +
  theme_classic() +
  #theme(aspect.ratio = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
        axis.title.y = element_text(size = 9, color = "black"),
        legend.position="none",
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(linewidth =0.2)) +
  labs(y = "Log2 Intensity") +
  xlab("")

ggsave(plot = Fig2B, width = 3, height = 3, dpi = 300, filename = "FIG2B.pdf")
