####################################################
# Analysis setup

# Path to the scripts
library("rstudioapi")

# the following line is for getting the path of your current open file
script_path <- getActiveDocumentContext()$path

dirpath <- dirname(script_path)

# The next line set the working directory to the relevant one:
setwd(dirname(dirpath))

# Load libraries
source("./scripts/libraries.R")
library("reshape2")

# load functions
source("./scripts/functions.R")

load("./RDS/tnbc_NoFib.RData")
dds_NoFib <- dds
rm(dds)

load("./RDS/deg_results_deseq2.RData")

####################################
# Data transformations and visualization
####################################
# Data transformations and visualization
norm_counts <- get.norm_dds(dds_NoFib)

norm_counts <- getHGNC(norm_counts)
colnames(norm_counts) <- gsub("TIL", "sTIL", colnames(norm_counts))

#################
# PCR Validation
#################
# Import PCR data
pcr <-
  read.csv(
    "./inputs/table.csv",
    sep = ",",
    header = TRUE,
    row.names = 1
  )

rownames(pcr) <- gsub("TIL", "sTIL", rownames(pcr))

pcr$Cells <-
  sapply(strsplit(as.character(row.names(pcr)), "-"), "[[", 2)


pcr <-
  pcr[, c(
    "CD28",
    "CCR7",
    "CD79B",
    "FCMR",
    "PAX5",
    "ELF3",
    "MAL2",
    "MUC1",
    "TFAP2A",
    "GPR37",
    "Cells"
  )]

library(reshape2)
library(ggpubr)

pcr.m <- melt(pcr, id.var = "Cells")

expt_plot <-
  ggplot(data = pcr.m, aes(x = variable, y = value)) + geom_boxplot(aes(fill =
                                                                          Cells)) +
  stat_compare_means(aes(group = Cells),
                     method = "t.test",
                     label = "p.format",
                     size = 2)

expt2 <-
  expt_plot + labs(x = "", y = "log2(normalized value) qRT-PCR") + scale_fill_manual(breaks = c("TC", "sTIL"),
                                                                                     values =
                                                                                       c("brown", "violet")) +
  theme(axis.text.x = element_text(
    face = "bold.italic",
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  )) +
  theme_classic() +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = 5.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "violet",
    alpha = 0.1
  )

#####################
# Boxplot Salmon
#####################

#Get ID of selected marker genes
tc_markers <-
  c("CD28",
    "CCR7",
    "CD79B",
    "FCMR",
    "PAX5",
    "ELF3",
    "MAL2",
    "MUC1",
    "TFAP2A",
    "GPR37")

#normalized data for selected ids
tc_markers_counts <- norm_counts[tc_markers, ]

tc_markers_counts <- t(tc_markers_counts)

tc_markers_counts <- as.data.frame(tc_markers_counts)

tc_markers_counts$Cells <-
  sapply(strsplit(as.character(row.names(tc_markers_counts)), "_"), "[[", 2)

tc_markers_counts.m <- melt(tc_markers_counts, id.var = "Cells")

rnaseq <-
  ggplot(data = tc_markers_counts.m, aes(x = variable, y = value)) + geom_boxplot(aes(fill =
                                                                                        Cells)) +
  stat_compare_means(aes(group = Cells),
                     method = "t.test",
                     label = "p.format",
                     size = 2)
salmon <-
  rnaseq + labs(x = "", y = "log2(normalized counts) RNA-seq") + scale_fill_manual(breaks = c("TC", "sTIL"),
                                                                                   values =
                                                                                     c("brown", "violet")) +
  theme(axis.text.x = element_text(
    face = "bold.italic",
    angle = 0,
    vjust = 0.5,
    hjust = 0.5
  )) +
  theme_classic() +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = 5.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "violet",
    alpha = 0.1
  )

#Multipanel figure 1
library(multipanelfigure)

figure2 <- multi_panel_figure(
  width = 200,
  height = 200,
  columns = 1,
  rows = 2,
  row_spacing = 5,
  column_spacing = 5
)

figure2 %<>% fill_panel(salmon, column = 1, row = 1)
figure2 %<>% fill_panel(expt2, column = 1, row = 2)

# Save figure
figure2 %>% save_multi_panel_figure(filename = "./images/Supplementary_Figure_6.tiff",
                                    dpi = 600,
                                    compression = "lzw")
