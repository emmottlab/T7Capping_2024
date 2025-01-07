# This script was created by Leandro Xavier Neves, PhD (https://orcid.org/0000-0002-6074-1025)

# Load packages
library("dplyr")
library("ggplot2")

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
int <- log2(int[,-1])

