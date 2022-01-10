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
# Data transformations (variance stabilization)
####################################################

# Excluding fibroblast
vsd_matrix_NoFib <- get.vst(dds_NoFib)
colnames(vsd_matrix_NoFib) <-
  gsub("TIL", "sTIL", colnames(vsd_matrix_NoFib))

####################################################
# Heatmap based on the Cheomokine signaling genes
####################################################

chemo <-
  read.table("./inputs/genesets_chemo.txt")

colnames(chemo) <- "chemokine"

# Get gene symbols
vsd_matrix_NoFib <- getHGNC(vsd_matrix_NoFib)

# Subset vsd matrix based on chemokine
deg_vsd_matrix <-
  subset(vsd_matrix_NoFib,
         rownames(vsd_matrix_NoFib) %in% chemo$chemokine)
dim(deg_vsd_matrix)

# Scale and center
deg_vsd_matrix_scaled <-
  t(scale(t(deg_vsd_matrix), center = TRUE, scale = TRUE))

# Reorder columns of a heatmap

colnames(deg_vsd_matrix_scaled)
deg_vsd_matrix_scaled <-
  deg_vsd_matrix_scaled [, c(
    "65_sTIL",
    "69_sTIL",
    "72_sTIL",
    "73_sTIL",
    "76_sTIL",
    "85_sTIL",
    "88_sTIL",
    "65_TC",
    "69_TC",
    "72_TC",
    "73_TC",
    "76_TC",
    "85_TC",
    "88_TC"
  )]

# Reorder rows of a heatmap
deg_vsd_matrix_scaled <-
  deg_vsd_matrix_scaled[chemo$chemokine,]

#Heatmap color

library(circlize)

mycols <-  colorRamp2(c(-1.0, 0, 1.0), c("blue", "white", "red"))

cm = ColorMapping(name = "",
                  col_fun = colorRamp2(
                    breaks = c(-1.0, 0, 1.0),
                    color = c("blue", "white", "red")
                  ))

grid.newpage()
color_mapping_legend(cm, title_gp = gpar(fontsize = 16))

# Annotation data frame

#reorder metadata
cluster_info <-
  meta_info_Coll[match(colnames(deg_vsd_matrix_scaled), rownames(meta_info_Coll)), ]

cluster_info <- as.data.frame(cluster_info)

annotation = data.frame(Cells = cluster_info$Cells)

ha = HeatmapAnnotation(
  df = annotation,
  col = list(Cells = c("TC" = "brown", "sTIL" = "violet")),
  na_col = "white",
  height = unit(0.25, "cm"),
  simple_anno_size_adjust = TRUE,
  show_legend = FALSE,
  annotation_legend_param = list(Cells = list(nrow = 1))
)

#Plot heatmap as grob
tiff(
  "./images/Figure_3B.tiff",
  width = 3,
  height = 8,
  units = 'in',
  res = 600,
  compression =  "lzw+p"
)
draw(
  Heatmap(
    deg_vsd_matrix_scaled,
    show_row_names = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_heatmap_legend = FALSE,
    col = mycols,
    column_names_gp = gpar(fontsize = 10),
    row_names_gp = gpar(fontsize = 10),
    top_annotation = ha
  )
)
dev.off()
