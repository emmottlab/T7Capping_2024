# Script for analysis of proteins derived from MNV subgenomic mRNA
# VP1, VP2, VF1.

#Read in log2 transformed data 
dat <- read.csv('mnv_T7_VP1_2DATA.csv')
head(dat)
View(dat)

#Load packages
install.packages("dplyr")
install.packages("ggplot")
install.packages("reshape")
install.packages("hrbrthemes")
install.packages("drawProteins")

library(dplyr)
library(ggplot2)
library(tidyverse)
library("dplyr")
library("ggplot2")
library("reshape")
library("hrbrthemes")
library("viridis")
library("stringr")
library("drawProteins")

#We want boxplot with two cell lines (CD300 and normal) with each having the Mock,WT and YGSN conditions plotted in one for VP1 


# Remove rows with NA in all samples
dat[,5:26][dat[,5:26] == 0] <- NA
dat <- dat[rowSums(is.na(dat[, 5:26])) != ncol(dat[, 5:26]), ]
View(dat)

# Plot box-plot
conditions2 <- gsub("_\\d+\\.Intensity","",colnames(dat[,5:26]))
view(dat)

long2 <- melt(dat)
long2$conditions <- gsub("_\\d+\\.Intensity","",long2$variable)
long2$conditions <- factor(long2$conditions, levels = unique(conditions2)) #converts the conditions column to a factor and specifies the desired order of the levels

# plot polyprotein first
long2_pp <- long2 %>% filter(grepl("Genome polyprotein", Description, fixed = TRUE))
long2_pp <- long2_pp %>% filter(!grepl("(", Description, fixed = TRUE))

Fig3B <- ggplot(long2_pp, aes(x = conditions2, y = value, fill = Description)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.2) +
  scale_fill_manual(values = rep(c("#666699"),13)) +
  scale_y_continuous(limits = c(15,27)) +
  scale_x_discrete(labels = c(gsub("BS_|CD_|BSCD300_", "",unique(long2_pp$conditions)))) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0) +
  theme_classic(base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 0, hjust = NULL, size = 6.5, color = "black"),
        axis.title.y = element_text(size = 7, color = "black"),
        legend.position="none",
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(linewidth =0.2)) +
  labs(y = expression("Log"[2] ~ "protein intensity)"))  +
  xlab("")

Fig3B

ggsave(plot = Fig3B, width = 3, height = 2.5, dpi = 300, filename = "FIG3B.svg")


#vp1
long3_pp <- long2 %>% filter(grepl("Capsid protein VP1", Description, fixed = TRUE))
long3_pp <- long3_pp %>% filter(!grepl("(", Description, fixed = TRUE))

conditions <- c(rep(paste0("Cond_0", 2:9), each = 5),rep(paste0("Cond_", 10:13), each = 5),rep("Mock", each = 5))
conditions <- conditions[c(61:65,1:60)]
group_names <- c(rep("Mock",5),rep(c(rep("WT", each = 5), rep("YGSN", each = 5)), 6))
  
FigVP1c2 <- ggplot(long3_pp, aes(x = conditions2, y = value, fill = Description)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.2) +
  scale_fill_manual(values = rep(c("#666699"),13)) +
  scale_y_continuous(limits = c(12,27.5)) +
  scale_x_discrete(labels = c(gsub("BS_|CD_|BSCD300_", "",unique(long3_pp$conditions)))) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0) +
  theme_classic(base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 0, hjust = NULL, size = 6.5, color = "black"),
        axis.title.y = element_text(size = 7, color = "black"),
        legend.position="none",
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(linewidth =0.2)) +
  labs(y = expression("Log"[2] ~ "protein intensity)"))  +
  xlab("")

FigVP1c2

ggsave(plot = FigVP1, width = 3, height = 2.5, dpi = 300, filename = "FIGVP1.svg")
  
  
#VP2
long4_pp <- long2 %>% filter(grepl("ORF3 Minor Capsid VP2", Description, fixed = TRUE))
long4_pp <- long4_pp %>% filter(!grepl("(", Description, fixed = TRUE))

conditions <- c(rep(paste0("Cond_0", 2:9), each = 5),rep(paste0("Cond_", 10:13), each = 5),rep("Mock", each = 5))
conditions <- conditions[c(61:65,1:60)]
group_names <- c(rep("Mock",5),rep(c(rep("WT", each = 5), rep("YGSN", each = 5)), 6))

FigVP2 <- ggplot(long4_pp, aes(x = conditions2, y = value, fill = Description)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.2) +
  scale_fill_manual(values = rep(c("#666699"),13)) +
  scale_y_continuous(limits = c(10,30)) +
  scale_x_discrete(labels = c(gsub("BS_|CD_|BSCD300_", "",unique(long4_pp$conditions)))) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0) +
  theme_classic(base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 0, hjust = NULL, size = 6.5, color = "black"),
        axis.title.y = element_text(size = 7, color = "black"),
        legend.position="none",
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(linewidth =0.2)) +
  labs(y = expression("Log"[2] ~ "protein intensity)"))  +
  xlab("")

FigVP2

ggsave(plot = FigVP1, width = 3, height = 2.5, dpi = 300, filename = "FIGVP2.svg")


#vIRULENCE FACTOR 1

long5_pp <- long2 %>% filter(grepl("Virulence factor 1", Description, fixed = TRUE))
long5_pp <- long5_pp %>% filter(!grepl("(", Description, fixed = TRUE))

conditions <- c(rep(paste0("Cond_0", 2:9), each = 5),rep(paste0("Cond_", 10:13), each = 5),rep("Mock", each = 5))
conditions <- conditions[c(61:65,1:60)]
group_names <- c(rep("Mock",5),rep(c(rep("WT", each = 5), rep("YGSN", each = 5)), 6))

FigVIRF1b <- ggplot(long5_pp, aes(x = conditions2, y = value, fill = Description)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.2) +
  scale_fill_manual(values = rep(c("#666699"),13)) +
  scale_y_continuous(limits = c(15,27.5)) +
  scale_x_discrete(labels = c(gsub("BS_|CD_|BSCD300_", "",unique(long5_pp$conditions)))) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0) +
  theme_classic(base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 0, hjust = NULL, size = 6.5, color = "black"),
        axis.title.y = element_text(size = 7, color = "black"),
        legend.position="none",
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(linewidth =0.2)) +
  labs(y = expression("Log"[2] ~ "(protein intensity)"))  +
  xlab("")

FigVIRF1b

ggsave(plot = FigVP1, width = 3, height = 2.5, dpi = 300, filename = "FIGVIRF1.svg")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
