####################################################
# Analysis setup
####################################################

# Path to the scripts
library("rstudioapi")

# the following line is for getting the path of your current open file
script_path <- getActiveDocumentContext()$path

dirpath <- dirname(script_path)

# The next line set the working directory to the relevant one:
setwd(dirname(dirpath))

# Load libraries
source("./scripts/libraries.R")
# load functions
source("./scripts/functions.R")

# Load DESEQ2 Data

# Without Fib for downstream analysis
load("./RDS/tnbc_NoFib.RData")

# Updated names
meta_info <- apply(meta_info, 2, function(y)
  gsub("TIL", "sTIL", y))
rownames(meta_info) <-
  gsub("TIL", "sTIL", rownames(meta_info))
meta_info <- as.data.frame(meta_info)


dds_NoFib <- dds
meta_info_NoFib <- meta_info
rm(dds)
rm(meta_info)

# Deseq2 results for heatmap
load("./RDS/deg_results_deseq2.RData")

# Updated names
meta_info_Coll <- apply(meta_info_Coll, 2, function(y)
  gsub("TIL", "sTIL", y))
rownames(meta_info_Coll) <-
  gsub("TIL", "sTIL", rownames(meta_info_Coll))
meta_info_Coll <- as.data.frame(meta_info_Coll)


####################################################
# GSVA analysis
####################################################

# Collapsed technical replicates
ddsColl <- collapseReplicates(dds_NoFib, dds_NoFib$Group)

# removing rows in which there are very few reads
# minimum 10 reads in 7 samples(7 because smallest group size)

keep <- rowSums(counts(ddsColl) >= 10) >= 7
table(keep)
dds_filter <- ddsColl[keep,]

topMatrix <- assay(dds_filter)
colnames(topMatrix) <- gsub("TIL", "sTIL", colnames(topMatrix))


#Convert the ensemble names to Entrez IDs (eliminate non-matches where necessary)
topMatrix_ent <- getEntrz(topMatrix)

#C2Broadset
c2 <- readRDS("./RDS/c2list.rds")

gene.sets <- c2$c2list.ez

# Truncate gene set names
names(gene.sets) <-
  ifelse(nchar(names(gene.sets)) > 40, paste0(strtrim(names(gene.sets), 35), '...'), names(gene.sets))

es.max <-
  gsva(
    topMatrix_ent,
    gene.sets,
    mx.diff = FALSE,
    verbose = FALSE,
    parallel.sz = 1
  )
es.dif <-
  gsva(
    topMatrix_ent,
    gene.sets,
    mx.diff = TRUE,
    verbose = FALSE,
    parallel.sz = 1
  )

par(mfrow = c(1, 2), mar = c(4, 4, 4, 1))

plot(
  density(as.vector(es.max)),
  main = "Maximum deviation from zero",
  xlab = "GSVA score",
  lwd = 2,
  las = 1,
  xaxt = "n",
  xlim = c(-0.75, 0.75),
  cex.axis = 0.8
)

axis(
  1,
  at = seq(-0.75, 0.75, by = 0.25),
  labels = seq(-0.75, 0.75, by = 0.25),
  cex.axis = 0.8
)

plot(
  density(as.vector(es.dif)),
  main = "Difference between largest\npositive and negative deviations",
  xlab = "GSVA score",
  lwd = 2,
  las = 1,
  xaxt = "n",
  xlim = c(-0.75, 0.75),
  cex.axis = 0.8
)

axis(
  1,
  at = seq(-0.75, 0.75, by = 0.25),
  labels = seq(-0.75, 0.75, by = 0.25),
  cex.axis = 0.8
)

p <- nrow(topMatrix_ent)    ## number of genes
n <- ncol(topMatrix_ent)    ## number of samples

nGS <- 100    ## number of gene sets
min.sz <- 10  ## minimum gene set size
max.sz <- 999 ## maximum gene set size

nGrp1 <- 7 ## number of samples in group 1
nGrp2 <- n - nGrp1 ## number of samples in group 2

#Perform GSVA
topMatrixGSVA <- gsva(
  topMatrix_ent,
  gene.sets,
  min.sz = min.sz,
  max.sz = max.sz,
  kcdf = "Poisson",
  abs.ranking = FALSE,
  verbose = TRUE
)

#Use limma to identify enriced pathways between TC and TILS
library(limma)

design <-
  model.matrix(~ factor(meta_info_Coll$Cells, levels = c("TC", "sTIL")))

colnames(design) <- c("Intercept", "sTILVsTC")

## fit linear model
fit <- lmFit(topMatrixGSVA, design)

## estimate moderated t-statistics
fit <- eBayes(fit)

sigPathways <-
  topTable(
    fit,
    coef = "sTILVsTC",
    number = Inf,
    p.value = 0.05,
    adjust = "BH"
  )
sigPathways <- sigPathways[abs(sigPathways$logFC) > 0.58,]

rownames(sigPathways) <- sigPathways$ID

# Write Sig pathways
#Create an object that can easily be written to disc
wObject <- data.frame(rownames(sigPathways), sigPathways)
colnames(wObject) <-
  c(
    "Pathway",
    "Log2FoldChange",
    "MeanExpression",
    "tStat",
    "Pvalue",
    "AdjustedPval",
    "Bvalue"
  )

#write.xlsx(
#  wObject,
#  row.names = TRUE,
#  file = paste("./GSVA/GSVA_sig_pathways", Sys.Date(), ".xlsx", sep = "_")
#)

#Filter the GSVA object to only include significant pathways
topMatrixGSVA_all <- topMatrixGSVA[rownames(head(sigPathways, 30)),]

heat <- t(scale(t(topMatrixGSVA_all)))

################
#GSVA heatmap
##################
#Heatmap color

library(circlize)

mycols <-  colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

cm = ColorMapping(name = "",
                  col_fun = colorRamp2(
                    breaks = c(-1.5, 0, 1.5),
                    color = c("blue", "white", "red")
                  ))

grid.newpage()
color_mapping_legend(cm, title_gp = gpar(fontsize = 16))

##################
# Annotation data frame

#reorder metadata
gsea_info <-
  meta_info_Coll[match(colnames(heat), rownames(meta_info_Coll)), ]

annotation = data.frame(Cells = gsea_info$Cells)

ha = HeatmapAnnotation(
  df = annotation,
  col = list(Cells = c("TC" = "brown", "sTIL" = "violet")),
  na_col = "white",
  height = unit(0.25, "cm"),
  simple_anno_size_adjust = TRUE,
  show_legend = FALSE,
  annotation_legend_param = list(Cells = list(nrow = 1))
)

# plot heatmap as grob object
#Plot heatmap as grob
tiff(
  "./images/Figure2_G.tiff",
  width = 8,
  height = 10,
  units = 'in',
  res = 600,
  compression =  "lzw+p"
)
draw(
  Heatmap(
    heat,
    #name = "x",
    show_row_dend = FALSE,
    width = unit(10, "cm"),
    height = unit(22, "cm"),
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    show_heatmap_legend = FALSE,
    col = mycols,
    row_title = "gene sets",
    top_annotation = ha
  )
  #annotation_legend_side = "bottom"
)
dev.off()
