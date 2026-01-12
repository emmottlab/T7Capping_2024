# This script was created by Leandro Xavier Neves, PhD (https://orcid.org/0000-0002-6074-1025)

# Load packages
library("dplyr")
library("ggplot2")
library("reshape")
library("hrbrthemes")
library("viridis")
library("stringr")
library("drawProteins")


# Set working directory
setwd("Y:/Leandro/2024/Projects/18_T7capping")

# Figure 2
# Load Fragpipe output files
data <- read.table("./01_timsTOF-HT_transfection_1/Fragpipe/LN_LFQ-MBR_trypsinKRnotP_1MC_skyline_noMBR/combined_protein.tsv",
                   sep = "\t", header = TRUE, quote = "")

# Subset MNV proteins
MNV <- dplyr::filter(data, Organism == "Norovirus (isolate Mouse/")

# Subset and log2 transform quantitative data "sample.Intensity" columns
quant_columns <- colnames(MNV[,210:274])
quant_columns <- sort(quant_columns)

int <- MNV[,c("Description",quant_columns)]
int[,quant_columns][int[,quant_columns] == 0] <- NA
int <- dplyr::mutate(int,across(where(is.numeric),log2))

# Plot box-plot
conditions <- c(rep(paste0("Cond_0", 2:9), each = 5),rep(paste0("Cond_", 10:13), each = 5),rep("Mock", each = 5))
conditions <- conditions[c(61:65,1:60)]
group_names <- c(rep("Mock",5),rep(c(rep("WT", each = 5), rep("YGSN", each = 5)), 6))

int_long <- melt(int)
int_long <- int_long[c(61:65,1:60),]
int_long <- cbind(int_long,conditions,group_names)
int_long$conditions <- factor(int_long$conditions, levels = unique(conditions)) #converts the conditions column to a factor and specifies the desired order of the levels

Fig2B <- ggplot(int_long, aes(x = conditions, y = value, fill = conditions)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.2) +
  stat_boxplot(geom = "errorbar", width = 0.3, coef = NULL, alpha = 0.8, linewidth = 0.2) +
  scale_fill_manual(values = rep(c("#666699"),13)) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  scale_y_continuous(limits = c(16,22)) +
  scale_x_discrete(labels = c("Mock",rep(c("WT","YGSN"), 6))) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0) +
  theme_classic(base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 0, hjust = NULL, size = 5, color = "black"),
        axis.title.y = element_text(size = 7, color = "black"),
        legend.position="none",
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(linewidth =0.2)) +
  labs(y = expression("Log"[2]~"protein intensity"))  +
  xlab("")

Fig2B

ggsave(plot = Fig2B, width = 3, height = 2.5, dpi = 300, filename = "FIG2B.png")

# Load Fragpipe combined_peptide output for MNV protein coverage map
peptide <- read.table("./01_timsTOF-HT_transfection_1/Fragpipe/LN_LFQ-MBR_trypsinKRnotP_1MC_skyline_noMBR/combined_peptide.tsv",
                      sep = "\t", header = TRUE, fill = TRUE, quote = "")

# Subset MNV peptides
peptide <- dplyr::filter(peptide, str_detect(peptide$Entry.Name, "_MNV") == TRUE)

#get UniProt features for MNV proteins identified
unique(peptide$Protein.ID) %>%
  drawProteins::get_features() %>%
  drawProteins::feature_to_dataframe() ->
  MNV_UniProt_data

#create canvas with protein length
p <- ggplot() +
  ylim(0.975, 1.3) +
  xlim(0,max(MNV_UniProt_data$end))

# add titles
p <- p + labs(x = "Amino acid number",
              y = "",)

#edit nsp names
POLG_ns <- MNV_UniProt_data[c(2,5:9),]
POLG_ns[1,2] <- "NS1/2 \n(p48)"  
POLG_ns[2,2] <- "NS3 \n(NTPase)"  
POLG_ns[3,2] <- "NS4 \n(p22)"  
POLG_ns[4,2] <- "NS5 \n(VPg)"  
POLG_ns[5,2] <- "NS6 \n(3CL)"  
POLG_ns[6,2] <- "NS7 \n(RdRP)"  

#plot chain with rounded corners
p <- p + statebins:::geom_rrect(data = POLG_ns[POLG_ns$type == "CHAIN",],
                                mapping=aes(xmin=begin,
                                            xmax=end,
                                            ymin=1-0.018,
                                            ymax=1+0.018),
                                radius = grid::unit(2, "pt"),
                                colour = "grey1", fill = "white")

#plot peptide coverage
p <- p + statebins:::geom_rrect(data = peptide,
                                mapping=aes(xmin=Start,
                                            xmax=End,
                                            ymin=1-0.0178,
                                            ymax=1+0.0178),
                                radius = grid::unit(2, "pt"),
                                alpha = 0.5)

#plot nsp names
p <- p + geom_label(data = POLG_ns[POLG_ns$type == "CHAIN",],
                    aes(x = begin + (end-begin)/2,
                        y = 1,
                        label = description),
                    fontface = "bold",
                    size = 3.8,
                    fill = NA,
                    label.size = NA)



# white background and remove y-axis
p <- p + theme_bw() +    # white background and text size
  theme() +
  #labs(fill="Domains") +
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank()) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  geom_rect(fill = NA) +
  theme(panel.border = element_blank(),axis.line.x = element_line(linetype = "solid", colour = "black")) +
  scale_x_continuous(breaks = seq(0,max(MNV_UniProt_data$end),150),limits = c(0,max(MNV_UniProt_data$end)))

p

ggsave(plot = p, dpi = 300, filename = "ORF1_cov.png")


# Figure 3
# Load Fragpipe output files
data2 <- read.table("./Deposit/LN_LFQ-MBR_trypsinKRnotP_1MC_skyline_noMBR_BSRT7complete_fiinal/combined_protein.tsv",
                   sep = "\t", header = TRUE, quote = "")

# Subset and normalise to total sum
quant <- data2[,81:102]

target_sum <- mean(colSums(quant))
quant_norm <- sweep(quant, 2, colSums(quant) / target_sum, FUN = "/")

# Add organism, protein description and gene to normalised data
quant_norm <- cbind(data2[,c("Protein.ID","Description","Gene","Organism")],quant_norm)

# Subset MNV proteins
MNV2 <- dplyr::filter(quant_norm, Organism == "Norovirus (isolate Mouse/")

# log2 transform
MNV2[,5:26][MNV2[,5:26] == 0] <- NA
MNV2 <- dplyr::mutate(MNV2,across(where(is.numeric),log2))

# Remove rows with NA in all samples
MNV2 <- MNV2[rowSums(is.na(MNV2[, 5:26])) != ncol(MNV2[, 5:26]), ]

# Plot box-plot
conditions2 <- gsub("_\\d+\\.Intensity","",colnames(MNV2[,5:26]))

long2 <- melt(MNV2)
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

ggsave(plot = Fig3B, width = 3, height = 2.5, dpi = 300, filename = "FIG3B.png")

# coverage map Figure 3
# Load Fragpipe combined_peptide output for MNV protein coverage map
peptide2 <- read.table("./Deposit/LN_LFQ-MBR_trypsinKRnotP_1MC_skyline_noMBR_BSRT7complete_fiinal/combined_peptide.tsv",
                      sep = "\t", header = TRUE, fill = TRUE, quote = "")

# Subset MNV peptides
peptide2 <- dplyr::filter(peptide2, str_detect(peptide2$Entry.Name, "POLG") == TRUE)

#get UniProt features for MNV polyprotein identified
c("Q80J95") %>%
  drawProteins::get_features() %>%
  drawProteins::feature_to_dataframe() ->
  MNV_UniProt_data

#create canvas with protein length
p2 <- ggplot() +
  ylim(0.975, 1.3) +
  xlim(0,max(MNV_UniProt_data$end))

# add titles
p2 <- p2 + labs(x = "Amino acid number",
              y = "",)

#plot chain with rounded corners
p2 <- p2 + statebins:::geom_rrect(data = POLG_ns[POLG_ns$type == "CHAIN",],
                                mapping=aes(xmin=begin,
                                            xmax=end,
                                            ymin=1-0.018,
                                            ymax=1+0.018),
                                radius = grid::unit(2, "pt"),
                                colour = "grey1", fill = "white")

#plot peptide coverage
p2 <- p2 + statebins:::geom_rrect(data = peptide2,
                                mapping=aes(xmin=Start,
                                            xmax=End,
                                            ymin=1-0.0178,
                                            ymax=1+0.0178),
                                radius = grid::unit(2, "pt"),
                                alpha = 0.5)

#plot nsp names
p2 <- p2 + geom_label(data = POLG_ns[POLG_ns$type == "CHAIN",],
                    aes(x = begin + (end-begin)/2,
                        y = 1,
                        label = description),
                    fontface = "bold",
                    size = 3.8,
                    fill = NA,
                    linewidth = NA)



# white background and remove y-axis
p2 <- p2 + theme_bw() +    # white background and text size
  theme() +
  #labs(fill="Domains") +
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank()) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  geom_rect(fill = NA) +
  theme(panel.border = element_blank(),axis.line.x = element_line(linetype = "solid", colour = "black")) +
  scale_x_continuous(breaks = seq(0,max(MNV_UniProt_data$end),150),limits = c(0,max(MNV_UniProt_data$end)))

p2 

ggsave(plot = p2, dpi = 300, filename = "ORF1_cov2.png")


# Save environment
save.image(paste0(Sys.Date(),"_BoxPlot_proteinCoverageMap.RData"))

