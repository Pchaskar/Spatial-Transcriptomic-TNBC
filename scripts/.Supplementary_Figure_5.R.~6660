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

# Load DESEQ2 Data
load("./RDS/tnbc_Fib.RData")

###
# Data transformations and visualization

vsd_matrix<-get.vst_rep(dds)


#PCA with Technical replicate L1 L2.
#Including fibroblast.
#Filtered.

# perform a PCA on the data in assay(x) for the selected genes.
pca <- prcomp(t(vsd_matrix))

# all PCs
df_out <- as.data.frame(pca$x)

# group based on cell type.

df_out$Cells <- sapply( strsplit(as.character(row.names(t(vsd_matrix))), "_"), "[[", 2 )# group on cell type TC, TIL, Fib.

df_out$lane <- sapply( strsplit(as.character(row.names(t(vsd_matrix))), "_"), "[[", 3 )# L1 and L2

head(df_out)

# percentage variance calculation.
percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=Cells, shape=lane))
p<-p+geom_point(size=3)+ xlab(percentage[1]) + ylab(percentage[2])
p<-p+geom_point(size=3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"))+
theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
      axis.title = element_text(size=10,face="bold"),
      legend.title = element_text(size=10,face="bold"),
      legend.text = element_text(size=10))

pca_ggplot<-p + scale_colour_manual(values = c("blue","brown","violet"))


tiff("./images/Supplementary_Figure_5_allgenes.tiff", width = 4, height = 3.5, 
     units = 'in', res = 600, compression =  "lzw+p")
pca_ggplot
dev.off()
