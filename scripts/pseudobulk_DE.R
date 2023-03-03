#######################################
# Load libraries
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(edgeR)
library(tidyverse)
library(magrittr)
library(purrr)
library(Matrix)
library(Matrix.utils)

###########################
# params
input_adata = 'data/PMCA_Jak2_raw_count.h5ad'
outdir = "data/"
test_condition = "Jak2_Homo"
ctrl_condition = "Jak2_WT"
outname_prefix = "PMCA_Jak2"

###########################
# Convert h5ad to h5seurat
input_h5seurat = str_replace(input_adata,".h5ad",".h5seurat")

if (!file.exists(input_h5seurat)) {
  Convert(input_adata, dest = "h5seurat")
}

#######################################
# Load data
seurat.obj = LoadH5Seurat(input_h5seurat)

# Normalize for later use
seurat.obj = NormalizeData(seurat.obj,
                           normalization.method = "RC",
                           scale.factor = 1e4)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = seurat.obj@assays$RNA@counts),
                            colData = seurat.obj@meta.data)

#######################################
# Aggregate counts per animal and celltype
groups <- colData(sce)[, c("celltype", "library")]
pb <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 

#######################################
# Create sample level metadata table
sids <- purrr::set_names(levels(sce$library))
n_cells <- as.numeric(table(sce$library))
m <- match(sids, sce$library)
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL)
ei = ei[,c("library","dataset","Condition","n_cells")]

#######################################
# Create pseudobulk metadata table
tab = table(sce$celltype, sce$library)
tab = as.data.frame(tab)
colnames(tab) = c("celltype","library","n_cells")
tab$pb_name = str_c(tab$celltype,"_",tab$library)
metadata <- left_join(tab, ei[, c("library", "dataset", "Condition")]) 

#######################################
# Split pseudo-bulk dataframe by cell type label
pb = as.data.frame(pb)
pb$pb_name = rownames(pb)
pb = left_join(pb, metadata[, c("pb_name", "celltype","library")]) 
pb = split(pb, pb$celltype)
pb = lapply(pb, function(x) {rownames(x) <- x$pb_name;
x$celltype <- NULL;
x$pb_name <- NULL;
x$library <- NULL; t(x) })

#######################################
# edgeR
celltypes <- levels(metadata$celltype)

for (celltype.to.comp in celltypes){
  comp = which(celltypes == celltype.to.comp)
  
  cluster_metadata <- metadata[which(metadata$celltype == celltypes[comp]), ]
  rownames(cluster_metadata) <- cluster_metadata$pb_name
  
  # Get pseudobulk counts
  counts <- pb[[celltypes[comp]]]
  
  pb_name_use = rownames(cluster_metadata)
  counts = counts[,pb_name_use]
  cluster_metadata = cluster_metadata[pb_name_use,]
  
  #######################################
  group = factor(cluster_metadata$Condition,levels = c(ctrl_condition, test_condition))
  batch = factor(cluster_metadata$dataset)
  
  #######################################
  # Filter low-expression genes
  seurat.obj.tmp = subset(seurat.obj, subset = celltype == celltype.to.comp)
  mean.CP10K = apply(GetAssayData(seurat.obj.tmp),1,mean)
  expressed = names(mean.CP10K)[which(mean.CP10K > 0.05)]
  counts.filt = counts[expressed,]
  
  #######################################
  y <- DGEList(counts=counts.filt,group=group)
  
  if (length(unique(batch))==1){
    design <- model.matrix(~group)
  } else if (min(table(group))==1){
    design <- model.matrix(~group)
  } else {
    design <- model.matrix(~batch+group)
  }
  design
  
  y <- calcNormFactors(y)
  y <- estimateDisp(y,design, robust=TRUE)

  fit <- glmFit(y,design)
  if (length(unique(batch))==1){
    lrt <- glmLRT(fit,coef=2)
  } else if (min(table(group))==1){
    lrt <- glmLRT(fit,coef=2)
  } else {
    n_batch = length(unique(batch))
    comp_coef = 1+n_batch
    lrt <- glmLRT(fit,coef=comp_coef)
  }
  
  write.table(topTags(lrt,n=nrow(lrt$table)),
              paste0(outdir,outname_prefix,"_",celltype.to.comp,"_DE_",test_condition,"_vs_",ctrl_condition,".txt"),
              row.names=TRUE,col.names=NA,sep="\t",quote=FALSE)
  
  res = as.data.frame(topTags(lrt,n=nrow(lrt$table)))
  resSig = res[abs(res$logFC)>0.5 & res$FDR<0.05,]
  write.table(resSig,
              paste0(outdir,outname_prefix,"_",celltype.to.comp,"_DE_",test_condition,"_vs_",ctrl_condition,"_Sig.txt"),
              row.names=TRUE,col.names=NA,sep="\t",quote=FALSE)
}

###########################
# Deleting intermediate files
file.remove(input_h5seurat)

cat("Analysis completed\n")
