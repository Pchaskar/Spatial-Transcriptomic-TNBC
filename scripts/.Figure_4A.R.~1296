####################################################
# Analysis setup

# Path to the scripts
library("rstudioapi")

# the following line is for getting the path of your current open file
script_path <- getActiveDocumentContext()$path 

dirpath<-dirname(script_path)

# The next line set the working directory to the relevant one:
setwd(dirname(dirpath))

# Load libraries
source("./scripts/libraries.R")
# load functions
source("./scripts/functions.R")

load("./RDS/tnbc_NoFib.RData")
dds_NoFib<-dds
rm(dds)

load("./RDS/deg_results_deseq2.RData")

####################################
# Data transformations and visualization

norm_counts<-get.norm_dds(dds_NoFib)

norm_counts<-getHGNC(norm_counts)

############################################
#Immnue cell proportion estimation
############################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)

# MCP counter on all samples

res_mcp_counter_all = deconvolute(norm_counts, "mcp_counter")

res_mcp_counter_all %>%
  gather(sample, score, -cell_type) %>%
  ggplot(aes(x=sample, y=score, color=cell_type)) +
  geom_point(size=4) +
  facet_wrap(~cell_type, scales="free_x", ncol=3) +
  scale_color_brewer(palette="Paired", guide=FALSE) +
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mcp_results<-res_mcp_counter_all %>% gather(sample, score, -cell_type)

mcp_results$sample<-as.character(mcp_results$sample)

mcp_results$ct<-mcp_results$sample

#rename tc and tils cell type
mcp_results$ct<-gsub(".*_", "\\1", mcp_results$ct) #65_TC-->TC

#sort mcp_results as ct 
mcp_results$sample <- factor(mcp_results$sample, levels = unique(mcp_results$sample[order(mcp_results$ct)]))

#Rename MCP counter output
colnames(mcp_results)<-c("cell_type","Samples","MCP-counter score","ct" )

#Remove monocytes
mcp_results_nomono<-mcp_results[! mcp_results$cell_type == 'Monocyte',]

#Rename cell types
immnue_cells<-c(`T cell`="T Cells", `B cell`="B-lineage", `T cell CD8+`="CD8+ T Cells", `Neutrophil`="Neutrophils",
                `NK cell`="NK", `Macrophage/Monocyte`="Monocytic lineage", `Myeloid dendritic cell`="Myeloid dendritic",
                `cytotoxicity score`="Cytotoxic lymphocytes", `Cancer associated fibroblast`="CAF", 
                `Endothelial cell`="Endothelial Cells")

mcp_results_nomono$immnue_cells <- as.character(immnue_cells[mcp_results_nomono$cell_type])

#Plot
tiff("./images/Figure4A.tiff", width = 7, height = 5, 
     units = 'in', res = 600, compression =  "lzw+p")
ggplot(mcp_results_nomono, aes(x=Samples, y=`MCP-counter score`, color=immnue_cells)) +
  geom_point(size=3) +
  facet_wrap(~immnue_cells, scales="free_x", ncol=5) +
  scale_color_brewer(palette="Paired", guide=FALSE) +
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6), 
        axis.text.y = element_text(size = 6),
        strip.text = element_text(size = 8)) 
dev.off()



