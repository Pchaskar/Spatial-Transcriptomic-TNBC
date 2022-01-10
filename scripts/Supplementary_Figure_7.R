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
library("reshape2")

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

###########################################
# Extra Figures and tables
#####################

#Get ID of selected marker genes
#co_stimulation_tcells markers

co_stimulation_markers<-c("CD2", "CD27", "CD28", "CD40LG", "CD226", 
                          "TNFSF14", "TNFRSF4", "TNFRSF8", "TNFRSF9","TNFRSF18", "TNFRSF25",
                           "ICOS","SLAMF1")

#normalized data for selected ids
co_stimulation_markers_counts<-norm_counts[co_stimulation_markers, ]

co_stimulation_markers_counts<-t(co_stimulation_markers_counts)

co_stimulation_markers_counts<-as.data.frame(co_stimulation_markers_counts)

co_stimulation_markers_counts$Cells<-sapply(strsplit(as.character(row.names(co_stimulation_markers_counts)), "_"), "[[", 2 )

co_stimulation_markers_counts.m <- melt(co_stimulation_markers_counts, id.var = "Cells")

co_stimulation_rnaseq<-ggplot(data = co_stimulation_markers_counts.m, aes(x=variable, y=value)) + geom_boxplot(aes(fill=Cells)) +
  stat_compare_means(aes(group = Cells), method = "t.test", label = "p.format", size= 2)
co_stimulation<-co_stimulation_rnaseq+labs(x ="", y = "log2(normalized counts) RNA-seq")+
  scale_fill_manual(breaks = c("TC", "TIL"), values=c("brown", "violet"))+
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold.italic", angle=40, vjust=0.5, hjust=0.5))
  
###################
#co_inhibition_T cells
###################

co_inhibition_markers<-c("BTLA", "CTLA4", "CD160", "CD244", "CD274", "HAVCR2",
                         "LAIR1","LAG3", "PDCD1LG2","TIGIT", "VSIR")     
                            
#normalized data for selected ids
co_inhibition_counts<-norm_counts[co_inhibition_markers, ]

co_inhibition_counts<-t(co_inhibition_counts)

co_inhibition_counts<-as.data.frame(co_inhibition_counts)

co_inhibition_counts$Cells<-sapply(strsplit(as.character(row.names(co_inhibition_counts)), "_"), "[[", 2 )

co_inhibition_counts.m <- melt(co_inhibition_counts, id.var = "Cells")

co_inhibition_rnaseq<-ggplot(data = co_inhibition_counts.m, aes(x=variable, y=value)) + geom_boxplot(aes(fill=Cells)) +
  stat_compare_means(aes(group = Cells), method = "t.test", label = "p.format", size= 2)
co_inhibition_Tcells<-co_inhibition_rnaseq+labs(x ="", y = "log2(normalized counts) RNA-seq")+
  scale_fill_manual(breaks = c("TC", "TIL"), values=c("brown", "violet"))+
  theme_classic()+
  theme(axis.text.x = element_text(face = "bold.italic", angle=40, vjust=0.5, hjust=0.5))
  

#Multipanel Figure
library(multipanelfigure)

figure3 <- multi_panel_figure(
  width = 200, height = 200,
  columns = 1, rows = 2, row_spacing=5,
  column_spacing=5)

figure3 %<>% fill_panel(co_stimulation, column = 1, row = 1)
figure3 %<>% fill_panel(co_inhibition_Tcells, column = 1, row = 2)


# Save figure
figure3 %>% save_multi_panel_figure(filename = "./images/Supplementary_Figure_6.tiff", dpi = 600, compression = "lzw")

