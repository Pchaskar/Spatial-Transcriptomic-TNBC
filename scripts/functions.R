# Functions
getHGNC<-function(data)
{
  geneid<-as.data.frame(apply(geneid,2,function(x)gsub('\\s+', '',x)))
  
  # add gene symbol
  data <- data[which(rownames(data) %in% rownames(geneid)),]
  geneid_data <- geneid[which(rownames(geneid) %in% rownames(data)),]
  geneid_data <- geneid_data[!duplicated(geneid_data$hgnc_symbol),]

  data <- data[match(rownames(geneid_data), rownames(data)),]
  rownames(data) <- geneid_data$hgnc_symbol
  return(data)
}

getHGNC_all<-function(data)
{
  require(biomaRt)
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  annots <- getBM(mart=mart,
                  attributes=c("ensembl_gene_id", "hgnc_symbol"),
                  filter="ensembl_gene_id",
                  values=rownames(data),
                  uniqueRows=TRUE)
  #annots <- annots[!duplicated(annots[,1]),]
  data <- data[which(rownames(data) %in% annots[,1]),]
  annots <- annots[which(annots[,1] %in% rownames(data)),]
  data <- data[match(annots[,1], rownames(data)),]
  data$symbol <- annots[,2]
  return(data)
}

getEntrz<-function(data)
{
  geneid<-as.data.frame(apply(geneid,2,function(x)gsub('\\s+', '',x)))
  
  # add gene symbol
  data <- data[which(rownames(data) %in% rownames(geneid)),]
  geneid_data <- geneid[which(rownames(geneid) %in% rownames(data)),]
  geneid_data <- geneid_data[!duplicated(geneid_data$entrezgene_id),]

  data <- data[match(rownames(geneid_data), rownames(data)),]
  rownames(data) <- geneid_data$entrezgene_id
  return(data)
}

##############
## Simple function to write .GMT (gene matrix transposed) format
## files.
## Get a list of genesets (ie list of character vectors)
## and write them to "filename"
write.gmt <- function(geneset.list, description.list=NA, filename)
{
  geneset.names <- names(geneset.list)
  
  if (is.na(description.list)) {
    write.set <- function(x) {
      write.table(t(c(geneset.names[[x]],"na",geneset.list[[x]])),
                  file=filename, col.names=FALSE, row.names=FALSE,
                  quote=FALSE, append=TRUE, sep="\t")
      cat("\n",file=filename,append=TRUE)
      paste("Geneset [",geneset.names[[x]],"] written to file", sep="")
    }
  } else {
    write.set <- function(x) {
      write.table(t(c(geneset.names[[x]],
                      description.list[[x]],
                      geneset.list[[x]])),
                  file=filename, col.names=FALSE, row.names=FALSE,
                  quote=FALSE, append=TRUE, sep="\t")
      cat("\n",file=filename,append=TRUE)
      paste("Geneset [",geneset.names[[x]],"] written to file", sep="")
    }
  }
  lapply(seq_along(geneset.list), write.set)
}

write.gct <- function(data, annot, filename)
{
  cat("#1.2\n",file=filename)
  cat(dim(data)[1:2],sep="\t",file=filename, append=TRUE)
  cat("\n",file=filename,append=TRUE)
  
  #outdf <- data.frame(t(data))
  outdf <- data.frame(data)
  
  outdf <- cbind("Description"=annot[dimnames(outdf)[[1]],]$Cells,outdf)
  outdf <- cbind("NAME"=rownames(outdf), outdf)
  rownames(outdf) <- NULL
  write.table(outdf,file=filename,append=TRUE,quote=FALSE,
              sep="\t",row.names=FALSE, na="")
}

call.gsea <- function(output.name="GSEA.analysis", iterations=1000)
{
  GSEA(
    ## IO
    input.ds = "aml.gct",
    input.cls = "aml.cls",
    gene.ann = "",
    gs.db = "aml.gmt",
    gs.ann = "",
    output.directory = "GSEA-output",
    doc.string = output.name,
    ## Base parameteres
    #non.interactive.run = F,
    reshuffling.type = "sample.labels",
    nperm = iterations,
    weighted.score.type = 1,
    nom.p.val.threshold = -1,
    fwer.p.val.threshold = -1,
    fdr.q.val.threshold = 0.25,
    topgs = 20,
    adjust.FDR.q.val = F,
    gs.size.threshold.min = 5,
    gs.size.threshold.max = 500,
    reverse.sign = F,
    preproc.type = 0,
    random.seed = 123456,
    ## Advanced parameters
    perm.type = 0,
    fraction = 1.0,
    replace = F,
    save.intermediate.results = F,
    #OLD.GSEA = F,
    use.fast.enrichment.routine = T)
}

get.cpm <- function(dds, meta_info)
{
  # Collapsed technical replicates
  ddsColl <- collapseReplicates(dds, dds$Group)
  meta_info_Coll<-meta_info[!duplicated(meta_info$Group), ]
  rownames(meta_info_Coll)<-meta_info_Coll$Group
  meta_info_Coll<-droplevels(meta_info_Coll)
  
  # Construct edgeR object
  
  se<-assay(ddsColl)
  
  coldata<-meta_info_Coll[match(colnames(se),rownames(meta_info_Coll)), ]
  
  library("edgeR")
  genetable <- data.frame(gene.id=rownames(se))
  y <- DGEList(counts=se, 
               samples=coldata, 
               genes=genetable)
  names(y)
  
  # DEsign formula same as deseq
  
  design <- model.matrix(~ Cells, y$samples)
  
  # Filter genes with low counts
  keep <- filterByExpr(y, design)
  table(keep)
  
  # recompute library siz
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  # normalization
  lcpm <- cpm(y, log = TRUE, prior.count = 0.25)
  
  topMatrix <- lcpm
  
  return(topMatrix)
}

get.cpm_rep <- function(dds, meta_info)
{
  
  # Construct edgeR object
  
  se<-assay(dds)
  
  coldata<-meta_info[match(colnames(se),rownames(meta_info)), ]
  
  library("edgeR")
  genetable <- data.frame(gene.id=rownames(se))
  y <- DGEList(counts=se, 
               samples=coldata, 
               genes=genetable)
  names(y)
  
  # DEsign formula same as deseq
  
  design <- model.matrix(~ Cells, y$samples)
  
  # Filter genes with low counts
  keep <- filterByExpr(y, design)
  table(keep)
  
  # recompute library siz
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  # normalization
  lcpm <- cpm(y, log = TRUE, prior.count = 0.25)
  
  topMatrix <- lcpm
  
  return(topMatrix)
}

get.vst <- function(dds)
{
  # Collapsed technical replicates
  ddsColl <- collapseReplicates(dds, dds$Group)
  
  
  # Pre-filtering the dataset
  # removing rows in which there are very few reads, we reduce the memory size of the dds data object
  # at least 4 samples with a count of 10 or higher
  
  keep <- rowSums(counts(ddsColl) >= 10) >= 4
  table(keep)
  dds_filter<-ddsColl[keep,]
  
  # VST
  vsd <- vst(dds_filter, blind=TRUE)
  vsd_matrix<-assay(vsd)
  
  return(vsd_matrix)
}

get.vst_rep <- function(dds)
{
  # Pre-filtering the dataset
  # removing rows in which there are very few reads, we reduce the memory size of the dds data object
  
  keep <- rowSums(counts(dds) >= 10) >= 4
  table(keep)
  dds_filter<-dds[keep,]
  
  # VST
  vsd <- vst(dds_filter, blind=TRUE)
  vsd_matrix<-assay(vsd)
  
  return(vsd_matrix)
}

get.norm_dds <- function(dds)
{
  # Collapsed technical replicates
  ddsColl <- collapseReplicates(dds, dds$Group)
  
  # Pre-filtering the dataset
  # removing rows in which there are very few reads, we reduce the memory size of the dds data object
  
  #keep <- rowSums(counts(ddsColl) >= 10) >= 7
  keep <- rowSums(counts(ddsColl) >= 10) >= 4
  
  table(keep)
  dds_filter<-ddsColl[keep,]
  
  # estimateSizeFactors
  dds_filter <- estimateSizeFactors(dds_filter)
  
  # normalized_counts
  normalized_counts <- counts(dds_filter, normalized=TRUE)
  
  normalized_counts<-log2(normalized_counts+1)
  
  return(normalized_counts)
}

get.tmm <- function(dds, meta_info)
{
  # Collapsed technical replicates
  ddsColl <- collapseReplicates(dds, dds$Group)
  meta_info_Coll<-meta_info[!duplicated(meta_info$Group), ]
  rownames(meta_info_Coll)<-meta_info_Coll$Group
  meta_info_Coll<-droplevels(meta_info_Coll)
  
  # Construct edgeR object
  
  se<-assay(ddsColl)
  
  coldata<-meta_info_Coll[match(colnames(se),rownames(meta_info_Coll)), ]
  
  library("edgeR")
  genetable <- data.frame(gene.id=rownames(se))
  y <- DGEList(counts=se, 
               samples=coldata, 
               genes=genetable)
  names(y)
  
  # DEsign formula same as deseq
  
  design <- model.matrix(~ Cells, y$samples)
  
  # Filter genes with low counts
  keep <- filterByExpr(y, design)
  table(keep)
  
  # recompute library siz
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  # normalization
  tmm <- calcNormFactors(y, method = "TMM")
  
  lcpm <- cpm(tmm, log = TRUE, prior.count = 0.25)
  
  topMatrix <- lcpm
  
  return(topMatrix)
}

get.tmm_mod_filter <- function(dds, meta_info)
{
  # Collapsed technical replicates
  ddsColl <- collapseReplicates(dds, dds$Group)
  meta_info_Coll<-meta_info[!duplicated(meta_info$Group), ]
  rownames(meta_info_Coll)<-meta_info_Coll$Group
  meta_info_Coll<-droplevels(meta_info_Coll)
  
  # Construct edgeR object
  
  se<-assay(ddsColl)
  
  coldata<-meta_info_Coll[match(colnames(se),rownames(meta_info_Coll)), ]
  
  library("edgeR")
  genetable <- data.frame(gene.id=rownames(se))
  y <- DGEList(counts=se, 
               samples=coldata, 
               genes=genetable)
  names(y)
  
  # DEsign formula same as deseq
  
  design <- model.matrix(~ Cells, y$samples)
  
  # Filter genes with low counts
  keep <- rowSums(cpm(y) > 0.5) >= 1
  table(keep)
  
  # recompute library siz
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  # normalization
  tmm <- calcNormFactors(y, method = "TMM")
  
  lcpm <- cpm(tmm, log = TRUE, prior.count = 0.25)
  
  topMatrix <- lcpm
  
  return(topMatrix)
}
