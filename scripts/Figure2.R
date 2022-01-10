####################################################
# Analysis setup
####################################################

# Path to the scripts
library("rstudioapi")
library("limma")

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

# With Fib for PCA
load("./RDS/tnbc_Fib.RData")
dds_Fib<-dds
meta_info_Fib<-meta_info
rm(dds)
rm(meta_info)

# Without Fib for downstream analysis
load("./RDS/tnbc_NoFib.RData")
dds_NoFib<-dds
meta_info_NoFib<-meta_info
rm(dds)
rm(meta_info)

# Deseq2 results for heatmap
load("./RDS/deg_results_deseq2.RData")

####################################################
# Data transformations (variance stabilization)
####################################################

# Including fibroblast
vsd_matrix_fib<-get.vst(dds_Fib)

# Excluding fibroblast
vsd_matrix_NoFib <- get.vst(dds_NoFib)

####################################################
#PCA with Technical replicate L1 L2 collapsed
####################################################

#Including fibroblast

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(vsd_matrix_fib))

# all PCs
df_out <- as.data.frame(pca$x)

# group based on cell type
df_out$Cells <- sapply( strsplit(as.character(row.names(t(vsd_matrix_fib))), "_"), "[[", 2 )# group on cell type TC, TIL, Fib
head(df_out)

# percentage variance calculation
percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=Cells))
p<-p+geom_point(size=3)+ xlab(percentage[1]) + ylab(percentage[2])
p<-p+geom_point(size=3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title = element_text(size=12,face="bold"),
        legend.title = element_text(size=12,face="bold"),
        legend.text = element_text(size=12))

pca_ggplot<-p + scale_colour_manual(values = c("blue","brown","violet"))

####################################################
# Heatmap based on the DEG
####################################################

sig_deg<-Sorted_res_diff[Sorted_res_diff$padj < 0.05, ]
sig_deg<-as.data.frame(sig_deg)

# Subset vsd matrix based on DEG 
deg_vsd_matrix <- subset(vsd_matrix_NoFib, rownames(vsd_matrix_NoFib) %in% rownames(sig_deg))
dim(deg_vsd_matrix)

deg_vsd_matrix<-getHGNC(deg_vsd_matrix)

# Scale and center
deg_vsd_matrix_scaled<-t(scale(t(deg_vsd_matrix), center = TRUE, scale = TRUE))

# Distance calculation based on correlation
c <- cor((deg_vsd_matrix_scaled), method="pearson")
d <- as.dist(1-c)

#Clustering of sample based on correlation distances
hc <- hclust(d, method = "complete", members=NULL)

#Distance calculation based on correlation
c_g <- cor(t(deg_vsd_matrix_scaled), method="pearson")
d_g <- as.dist(1-c_g)

#Clustering of sample based on correlation distances
hr <- hclust(d_g, method = "complete", members=NULL)


#Heatmap color

library(circlize)

mycols <-  colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

cm = ColorMapping(name = "",
                  col_fun = colorRamp2(breaks=c(-1.5, 0, 1.5), color=c("blue", "white", "red")))

grid.newpage()
color_mapping_legend(cm, title_gp = gpar(fontsize = 16))

# Annotation data frame

#reorder metadata
cluster_info<-meta_info_Coll[match(colnames(deg_vsd_matrix_scaled),rownames(meta_info_Coll)), ]

annotation = data.frame(Cells=cluster_info$Cells)

ha = HeatmapAnnotation(df = annotation,
                       col = list(Cells=c("TC" = "brown", "TIL" = "violet")),
                       na_col = "white",
                       height = unit(0.25, "cm"),
                       simple_anno_size_adjust = TRUE,
                       show_legend = FALSE,
                       annotation_legend_param = list(
                         Cells = list(nrow = 1))
)

# Top DEG labels
sig_deg<-getHGNC(sig_deg)

sig_deg_up<-sig_deg[sig_deg$log2FoldChange > 0, ]
sig_deg_up_20top<-head(sig_deg_up[order(sig_deg_up$padj, decreasing=FALSE), ], 20)

sig_deg_down<-sig_deg[sig_deg$log2FoldChange < 0, ]
sig_deg_down_20top<-head(sig_deg_down[order(sig_deg_down$padj, decreasing=FALSE), ], 20)

sig_top20<- rbind(sig_deg_up_20top, sig_deg_down_20top)

at_up = which(rownames(deg_vsd_matrix_scaled) %in% rownames(sig_deg_up_20top))
at_down = which(rownames(deg_vsd_matrix_scaled) %in% rownames(sig_deg_down_20top))

at_sig=which(rownames(deg_vsd_matrix_scaled) %in% rownames(sig_top20))

ha_de_up = rowAnnotation(link = anno_mark(at_up, 
                                          labels = rownames(deg_vsd_matrix_scaled[at_up,]), 
                                          labels_gp = gpar(fontsize = 5, col="red")
                                          ))

ha_de_down = rowAnnotation(link = anno_mark(at_down, 
                                            labels = rownames(deg_vsd_matrix_scaled[at_down,]), 
                                            labels_gp = gpar(fontsize = 5, col="blue")
                                            ))

ha_top40 = rowAnnotation(link = anno_mark(at_sig, 
                                            labels = rownames(deg_vsd_matrix_scaled[at_sig,]), 
                                            labels_gp = gpar(fontsize = 5, 
            col=c(rep("red", 3), rep("blue", 7) ,rep("red", 7), rep("blue", 5), rep("red", 3),
            rep("blue", 4), rep("red", 1), rep("blue", 1), rep("red", 1),
            rep("blue", 1), rep("red", 3), rep("blue", 2), rep("red", 2))
)))

#Plot heatmap as grob
exp_heatmap<-grid::grid.grabExpr(draw(Heatmap(deg_vsd_matrix_scaled,
                                              show_row_names = FALSE, show_row_dend = FALSE, 
                                              cluster_rows = hr,
                                              cluster_columns = hc,
                                              heatmap_legend_param = list(title = "scale", 
                                                                          color_bar = "continuous"),
                                              col=mycols, column_names_gp = gpar(fontsize = 7), row_title="genes",
                                              top_annotation =ha)+ha_top40))

####################################################
# GSVA analysis
####################################################

# Collapsed technical replicates
ddsColl <- collapseReplicates(dds_NoFib, dds_NoFib$Group)

# removing rows in which there are very few reads
# minimum 10 reads in 7 samples(7 because smallest group size)

keep <- rowSums(counts(ddsColl) >= 10) >= 7
table(keep)
dds_filter<-ddsColl[keep,]

topMatrix <- assay(dds_filter)

#Convert the ensemble names to Entrez IDs (eliminate non-matches where necessary)
topMatrix_ent<-getEntrz(topMatrix)

#C2Broadset
c2<-readRDS("./RDS/c2list.rds")

gene.sets <- c2$c2list.ez

# Truncate gene set names
names(gene.sets)<-ifelse(nchar(names(gene.sets)) > 40, paste0(strtrim(names(gene.sets), 35), '...'), names(gene.sets))

es.max <- gsva(topMatrix_ent, gene.sets, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
es.dif <- gsva(topMatrix_ent, gene.sets, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)

par(mfrow=c(1,2), mar=c(4, 4, 4, 1))

plot(density(as.vector(es.max)), main="Maximum deviation from zero",
     xlab="GSVA score", lwd=2, las=1, xaxt="n", xlim=c(-0.75, 0.75), cex.axis=0.8)

axis(1, at=seq(-0.75, 0.75, by=0.25), labels=seq(-0.75, 0.75, by=0.25), cex.axis=0.8)

plot(density(as.vector(es.dif)), main="Difference between largest\npositive and negative deviations",
     xlab="GSVA score", lwd=2, las=1, xaxt="n", xlim=c(-0.75, 0.75), cex.axis=0.8)

axis(1, at=seq(-0.75, 0.75, by=0.25), labels=seq(-0.75, 0.75, by=0.25), cex.axis=0.8)

p <- nrow(topMatrix_ent)    ## number of genes
n <- ncol(topMatrix_ent)    ## number of samples

nGS <- 100    ## number of gene sets
min.sz <- 10  ## minimum gene set size
max.sz <- 999 ## maximum gene set size

nGrp1 <- 7 ## number of samples in group 1
nGrp2 <- n - nGrp1 ## number of samples in group 2

#Perform GSVA   
topMatrixGSVA <- gsva(topMatrix_ent, 
                      gene.sets, 
                      min.sz=min.sz, 
                      max.sz=max.sz,
                      kcdf="Poisson",
                      abs.ranking=FALSE, 
                      verbose=TRUE)

#Use limma to identify enriced pathways between TC and TILS

design <- model.matrix(~ factor(meta_info_Coll$Cells, levels=c("TC", "TIL")))

colnames(design) <- c("Intercept", "TILVsTC")

## fit linear model
fit <- lmFit(topMatrixGSVA, design)

## estimate moderated t-statistics
fit <- eBayes(fit)

sigPathways <- topTable(fit, coef="TILVsTC", number=Inf, p.value=0.05, adjust="BH")
sigPathways <- sigPathways[abs(sigPathways$logFC)>0.58,]

rownames(sigPathways)<-sigPathways$ID

#Pathways TC, logFCcutoff <- log2(1.5)
#sigPathways_TC <- sigPathways[sigPathways$logFC < 0.58,]

#Pathway TILS, logFCcutoff <- log2(1.5)
#sigPathways_TIL <- sigPathways[sigPathways$logFC > 0.58,]

# Write Sig pathways
#Create an object that can easily be written to disc
wObject <- data.frame(rownames(sigPathways), sigPathways)
colnames(wObject) <- c("Pathway","Log2FoldChange","MeanExpression","tStat","Pvalue","AdjustedPval","Bvalue")

write.xlsx(wObject,row.names = TRUE, file = paste("./GSVA/GSVA_sig_pathways",Sys.Date(), ".xlsx", sep="_"))

#Filter the GSVA object to only include significant pathways
topMatrixGSVA_all <- topMatrixGSVA[rownames(sigPathways),]

heat <- t(scale(t(topMatrixGSVA_all)))

################
#GSVA heatmap
##################
#Heatmap color

library(circlize)

mycols <-  colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

cm = ColorMapping(name = "",
                  col_fun = colorRamp2(breaks=c(-1.5, 0, 1.5), color=c("blue", "white", "red")))

grid.newpage()
color_mapping_legend(cm, title_gp = gpar(fontsize = 16))

##################
# Annotation data frame

#reorder metadata
gsea_info<-meta_info_Coll[match(colnames(heat),rownames(meta_info_Coll)), ]

annotation = data.frame(Cells=gsea_info$Cells)

ha = HeatmapAnnotation(df = annotation,
                       col = list(Cells=c("TC" = "brown", "TIL" = "violet")),
                       na_col = "white",
                       height = unit(0.25, "cm"),
                       simple_anno_size_adjust = TRUE,
                       show_legend = FALSE,
                       annotation_legend_param = list(
                         Cells = list(nrow = 1))
)

# plot heatmap as grob object
gsva_heatmap<-grid::grid.grabExpr(draw(Heatmap(heat,
                                               #name = "x",
                                               show_row_dend = FALSE,
                                               width = unit(3.5, "cm"), height = unit(16.5, "cm"),
                                               show_row_names = TRUE,
                                               row_names_gp = gpar(fontsize = 3),
                                               column_names_gp = gpar(fontsize = 7),
                                               show_heatmap_legend = FALSE,
                                               col=mycols,
                                               row_title="gene sets",
                                               top_annotation =ha)
                                       #annotation_legend_side = "bottom"
                                       ))

################
#Multipanelfigure
################

library(multipanelfigure)

# Layout
figure1 <- multi_panel_figure(
  width = 200, height = 200,
  columns = 6, rows = 7, row_spacing=4,
  column_spacing=3)


# Plot Multipanelfigure
figure1 %<>% fill_panel(pca_ggplot, column = 1:3, row = 1:2)
figure1 %<>% fill_panel(exp_heatmap, column = 1:3, row = 3:7, scaling="shrink")
figure1 %<>% fill_panel(gsva_heatmap, column = 4:6, row = 1:7, scaling="shrink")


# Save figure
figure1 %>% save_multi_panel_figure(filename = "./images/Figure3.tiff", dpi = 600, compression = "lzw")

